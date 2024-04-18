#ifndef CASM_clexmonte_IncrementalConditionsStateGenerator
#define CASM_clexmonte_IncrementalConditionsStateGenerator

#include <map>
#include <string>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/run/ConfigGenerator.hh"
#include "casm/clexmonte/run/RunData.hh"
#include "casm/clexmonte/run/StateGenerator.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/run_management/State.hh"

// io - completed_runs.json
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/misc/parse_array.hh"
#include "casm/clexmonte/run/io/json/RunData_json_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief Generates a series of states by constant conditions increments
///
/// The run information needed to check completion and generate subsequent
/// states is the vector of data from the previous runs.
///
/// This method generates states using the following steps:
/// 1) Set indepedently determined conditions, using:
///    \code
///    monte::ValueMap conditions = make_incremented_values(
///        initial_conditions, conditions_increment, completed_runs.size());
///    \endcode
/// 2) Generate an initial configuration, using:
///    \code
///    config_type configuration =
///        (dependent_runs && completed_runs.size()) ?
///            ? completed_runs.back().final_state->configuration
///            : config_generator(conditions, completed_runs);
///    \endcode
/// 3) Make a state, using:
///    \code
///    state_type state(configuration, conditions)
///    \endcode
/// 4) Apply custom state modifiers, using:
///    \code
///    for (auto const &f : m_modifiers) {
///      f(state);
///    }
///    \endcode
/// 5) Return `state`.
class IncrementalConditionsStateGenerator : public StateGenerator {
 public:
  /// \brief Constructor
  ///
  /// \param _config_generator Function to generate configurations for the
  ///     initial state from the indepedently determined conditions and the
  ///     data from previous runs.
  /// \param _initial_conditions The "indepedently determined conditions" for
  ///     the initial state.
  /// \param _conditions_increment The conditions to be changed between states,
  ///     and the amount to change them. A key in `_conditions_increment`
  ///     must also be a key of `_initial_conditions`.
  /// \param _n_states The total number of states to generate. Includes the
  ///     initial state.
  /// \param _dependent_runs If true, use the last configuration as the starting
  ///     point for the next state. If false, always use the configuration of
  ///     the initial state.
  /// \param _modifiers Functions that modify the generated state,
  ///     for instance to set the composition condition for canonical
  ///     calculations based on the composition of the generated or input
  ///     configuration so that it doesn't have to be pre-determined by
  ///     the user.
  IncrementalConditionsStateGenerator(
      std::shared_ptr<system_type> system, RunDataOutputParams output_params,
      std::unique_ptr<ConfigGenerator> _config_generator,
      monte::ValueMap const &_initial_conditions,
      monte::ValueMap const &_conditions_increment, Index _n_states,
      bool _dependent_runs,
      std::vector<StateModifyingFunction> const &_modifiers = {})
      : m_system(system),
        m_output_params(output_params),
        m_config_generator(std::move(_config_generator)),
        m_initial_conditions(_initial_conditions),
        m_conditions_increment(_conditions_increment),
        m_n_states(_n_states),
        m_dependent_runs(_dependent_runs),
        m_modifiers(_modifiers) {
    std::stringstream msg;
    msg << "Error constructing IncrementalConditionsStateGenerator: "
        << "Mismatch between initial conditions and conditions increment.";
    if (is_mismatched(m_initial_conditions, m_conditions_increment)) {
      throw std::runtime_error(msg.str());
    }
  }

  /// \brief Check if all requested runs have been completed
  bool is_complete() override { return m_completed_runs.size() >= m_n_states; }

  /// \brief Return the next state
  state_type next_state() override {
    if (m_dependent_runs && m_completed_runs.size() &&
        !m_completed_runs.back().final_state.has_value()) {
      throw std::runtime_error(
          "Error in IncrementalConditionsStateGenerator: when "
          "dependent_runs==true, must save the final state of the last "
          "completed run");
    }

    // Make conditions
    monte::ValueMap conditions = make_incremented_values(
        m_initial_conditions, m_conditions_increment, m_completed_runs.size());

    // Make configuration
    config_type configuration =
        (m_dependent_runs && m_completed_runs.size())
            ? m_completed_runs.back().final_state->configuration
            : (*m_config_generator)(conditions, m_completed_runs);

    // Make state
    state_type state(configuration, conditions);

    // Apply custom modifiers
    for (auto const &f : m_modifiers) {
      f(state, static_cast<monte::OccLocation *>(nullptr));
    }

    // Finished
    return state;
  }

  void push_back(RunData const &run_data) override {
    if (m_completed_runs.size() && !m_output_params.do_save_all_final_states) {
      m_completed_runs.back().final_state.reset();
    }
    m_completed_runs.push_back(run_data);
    if (!m_output_params.do_save_all_initial_states) {
      m_completed_runs.back().initial_state.reset();
    }
    if (!m_output_params.do_save_last_final_state &&
        !m_output_params.do_save_all_final_states) {
      m_completed_runs.back().final_state.reset();
    }
  }

  std::vector<RunData> const &completed_runs() const override {
    return m_completed_runs;
  }

  Index n_completed_runs() const override { return m_completed_runs.size(); }

  void read_completed_runs() override {
    m_completed_runs.clear();
    if (m_output_params.output_dir.empty()) {
      return;
    }

    fs::path completed_runs_path =
        m_output_params.output_dir / "completed_runs.json";
    if (!fs::exists(completed_runs_path)) {
      return;
    }

    jsonParser json(completed_runs_path);
    ParentInputParser parser{json};
    auto subparser = parser.parse_as_with<std::vector<RunData>>(
        parse_array<RunData, config::SupercellSet &>, *m_system->supercells);

    std::stringstream ss;
    ss << "Error in IncrementalConditionsStateGenerator: failed to read "
       << completed_runs_path;
    std::runtime_error error_if_invalid{ss.str()};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    m_completed_runs = *subparser->value;
  }

  void write_completed_runs() const override {
    if (m_output_params.output_dir.empty()) {
      return;
    }

    fs::path completed_runs_path =
        m_output_params.output_dir / "completed_runs.json";
    fs::create_directories(m_output_params.output_dir);
    SafeOfstream file;
    file.open(completed_runs_path);
    jsonParser json;
    to_json(m_completed_runs, json, m_output_params.write_initial_states,
            m_output_params.write_final_states);
    json.print(file.ofstream(), -1);
    file.close();
  }

 private:
  /// System data
  std::shared_ptr<system_type> m_system;

  RunDataOutputParams m_output_params;
  std::vector<RunData> m_completed_runs;
  std::unique_ptr<ConfigGenerator> m_config_generator;
  monte::ValueMap m_initial_conditions;
  monte::ValueMap m_conditions_increment;
  Index m_n_states;
  bool m_dependent_runs;
  std::vector<StateModifyingFunction> m_modifiers;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
