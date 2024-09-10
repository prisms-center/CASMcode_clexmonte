#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/io/json/EventState_json_io.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/clexmonte/methods/kinetic_monte_carlo.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/modifying_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

/// \brief Data calculated for a single event in a single state
struct EventState {
  bool is_allowed;      ///< Is allowed given current configuration
  bool is_normal;       ///< Is "normal" (dEa > 0.0) && (dEa > dEf)
  double dE_final;      ///< Final state energy, relative to initial state
  double Ekra;          ///< KRA energy
  double dE_activated;  ///< Activation energy, relative to initial state
  double freq;          ///< Attempt frequency
  double rate;          ///< Occurance rate
};

jsonParser &to_json(EventState const &event_state, jsonParser &json) {
  json["is_allowed"] = event_state.is_allowed;
  if (event_state.is_allowed) {
    json["is_normal"] = event_state.is_normal;
    json["dE_final"] = event_state.dE_final;
    json["Ekra"] = event_state.Ekra;
    json["dE_activated"] = event_state.dE_activated;
    json["freq"] = event_state.freq;
    json["rate"] = event_state.rate;
  }
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    PrimEventData const &prim_event_data) {
  to_json(prim_event_data, json);
  to_json(event_state, json);
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    EventData const &event_data,
                    PrimEventData const &prim_event_data) {
  json["unitcell_index"] = event_data.unitcell_index;
  json["linear_site_index"] = event_data.event.linear_site_index;
  to_json(event_state, json, prim_event_data);
  return json;
}

void print(std::ostream &out, EventState const &event_state) {
  out << "is_allowed: " << std::boolalpha << event_state.is_allowed
      << std::endl;
  if (event_state.is_allowed) {
    out << "dE_activated: " << event_state.dE_activated << std::endl;
    out << "dE_final: " << event_state.dE_final << std::endl;
    out << "is_normal: " << std::boolalpha << event_state.is_normal
        << std::endl;
    out << "Ekra: " << event_state.Ekra << std::endl;
    out << "freq: " << event_state.freq << std::endl;
    out << "rate: " << event_state.rate << std::endl;
  }
}

void print(std::ostream &out, EventState const &event_state,
           PrimEventData const &prim_event_data) {
  out << "prim_event_index: " << prim_event_data.prim_event_index << std::endl;
  out << "event_type_name: " << prim_event_data.event_type_name << std::endl;
  out << "equivalent_index: " << prim_event_data.equivalent_index << std::endl;
  out << "is_forward: " << std::boolalpha << prim_event_data.is_forward
      << std::endl;
  out << "occ_init: " << prim_event_data.occ_init << std::endl;
  out << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(out, event_state);
}

void print(std::ostream &out, EventState const &event_state,
           EventData const &event_data, PrimEventData const &prim_event_data) {
  out << "prim_event_index: " << prim_event_data.prim_event_index << std::endl;
  out << "unitcell_index: " << event_data.unitcell_index << std::endl;
  out << "event_type_name: " << prim_event_data.event_type_name << std::endl;
  out << "equivalent_index: " << prim_event_data.equivalent_index << std::endl;
  out << "is_forward: " << std::boolalpha << prim_event_data.is_forward
      << std::endl;
  out << "linear_site_index: " << event_data.event.linear_site_index
      << std::endl;
  out << "occ_init: " << prim_event_data.occ_init << std::endl;
  out << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(out, event_state);
}

/// \brief Event rate calculation for a particular KMC event
///
/// EventStateCalculator is used to separate the event calculation from the
/// event definition data in PrimEventData. All symmetrically equivalent
/// events can use the same EventStateCalculator, but a simple approach
/// is to create one for each distinct event associated with the primitive
/// cell.
class EventStateCalculator {
 public:
  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name)
      : m_system(_system), m_event_type_name(_event_type_name) {}

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state) {
    // supercell-specific
    m_state = state;
    if (m_state == nullptr) {
      throw std::runtime_error(
          "Error setting EventStateCalculator state: state is empty");
    }
    m_temperature = &m_state->conditions.scalar_values.at("temperature");
    m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

    // set and validate event clex
    LocalMultiClexData event_local_multiclex_data =
        get_local_multiclex_data(*m_system, m_event_type_name);
    m_event_clex = get_local_multiclex(*m_system, *m_state, m_event_type_name);
    std::map<std::string, Index> _glossary =
        event_local_multiclex_data.coefficients_glossary;

    auto _check_coeffs = [&](Index &coeff_index, std::string key) {
      if (!_glossary.count(key)) {
        std::stringstream ss;
        ss << "Error constructing " << m_event_type_name
           << " EventStateCalculator: No " << key << " cluster expansion";
        throw std::runtime_error(ss.str());
      }
      coeff_index = _glossary.at(key);
      if (coeff_index < 0 ||
          coeff_index >= m_event_clex->coefficients().size()) {
        std::stringstream ss;
        ss << "Error constructing " << m_event_type_name
           << " EventStateCalculator: " << key << " index out of range";
        throw std::runtime_error(ss.str());
      }
    };
    _check_coeffs(m_kra_index, "kra");
    _check_coeffs(m_freq_index, "freq");
  }

  /// \brief Pointer to current state
  state_type const *state() const { return m_state; }

  /// \brief Current state's reciprocal temperature
  double beta() const { return 1.0 / (CASM::KB * *this->m_temperature); }

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, EventData const &event_data,
                             PrimEventData const &prim_event_data) const {
    clexulator::ConfigDoFValues const *dof_values =
        m_formation_energy_clex->get();

    // Check if event is allowed based on the current occupation
    int i = 0;
    for (Index l : event_data.event.linear_site_index) {
      if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
        state.is_allowed = false;
        state.rate = 0.0;
        return;
      }
      ++i;
    }
    state.is_allowed = true;

    // calculate change in energy to final state
    state.dE_final = m_formation_energy_clex->occ_delta_value(
        event_data.event.linear_site_index, prim_event_data.occ_final);

    // calculate KRA and attempt frequency
    Eigen::VectorXd const &event_values = m_event_clex->values(
        event_data.unitcell_index, prim_event_data.equivalent_index);
    state.Ekra = event_values[m_kra_index];
    state.freq = event_values[m_freq_index];

    // calculate energy in activated state, check if "normal", calculate rate
    state.dE_activated = state.dE_final * 0.5 + state.Ekra;
    state.is_normal =
        (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
    if (state.dE_activated < state.dE_final)
      state.dE_activated = state.dE_final;
    if (state.dE_activated < 0.0) state.dE_activated = 0.0;
    state.rate = state.freq * exp(-this->beta() * state.dE_activated);
  }

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// State to use
  state_type const *m_state;

  /// Current state's temperature
  double const *m_temperature;

  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
  Index m_kra_index;
  Index m_freq_index;
};

/// \brief CompleteEventCalculator is an event calculator with the required
/// interface for the
///     classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct CompleteEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Prim event calculators - order must match prim_event_list
  std::vector<EventStateCalculator> const &prim_event_calculators;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Write to warn about non-normal events
  Log &event_log;

  // Note: to keep all event state calculations, comment out this:
  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Count not-normal events
  Index not_normal_count;

  CompleteEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<EventStateCalculator> const &_prim_event_calculators,
      std::map<EventID, EventData> const &_event_list,
      Log &_event_log = CASM::err_log())
      : prim_event_list(_prim_event_list),
        prim_event_calculators(_prim_event_calculators),
        event_list(_event_list),
        event_log(_event_log),
        not_normal_count(0) {}

  /// \brief Get CASM::monte::OccEvent corresponding to given event ID
  double calculate_rate(EventID const &id) {
    EventData const &event_data = event_list.at(id);
    PrimEventData const &prim_event_data =
        prim_event_list.at(id.prim_event_index);
    // Note: to keep all event state calculations, uncomment this:
    // EventState &event_state = event_data.event_state;
    prim_event_calculators.at(id.prim_event_index)
        .calculate_event_state(event_state, event_data, prim_event_data);

    // ---
    // can check event state and handle non-normal event states here
    // ---
    if (event_state.is_allowed && !event_state.is_normal) {
      event_log << "---" << std::endl;
      print(event_log.ostream(), event_state, event_data, prim_event_data);
      event_log << std::endl;
      ++not_normal_count;
    }

    return event_state.rate;
  }
};

/// \brief Data for kinetic Monte Carlo events
///
/// Includes:
/// - prim event list
/// - prim impact info
/// - event state calculators: one per prim event, given a state pointer and
///   can then calculate event energies, attempt frequency, and rate for the
///   for the current state on request
/// - complete event list
/// - CompleteEventCalculator: uses event state calculator and complete event
///   list to calculate a rate given an event ID
struct KineticEventData {
  KineticEventData(std::shared_ptr<system_type> _system) : system(_system) {
    if (!is_clex_data(*system, "formation_energy")) {
      throw std::runtime_error(
          "Error constructing KineticEventData: no 'formation_energy' clex.");
    }

    prim_event_list = clexmonte::make_prim_event_list(*system);
    prim_impact_info_list = clexmonte::make_prim_impact_info_list(
        *system, prim_event_list, {"formation_energy"});
  }

  /// \brief Update for given state, conditions, and occupants
  void update(state_type const &state, monte::OccLocation const &occ_location,
              std::vector<EventFilterGroup> const &event_filters) {
    // These are constructed/re-constructed so cluster expansions point
    // at the current state
    prim_event_calculators.clear();
    for (auto const &prim_event_data : prim_event_list) {
      prim_event_calculators.emplace_back(system,
                                          prim_event_data.event_type_name);
      prim_event_calculators.back().set(&state);
    }

    // TODO: rejection-clexmonte option does not require impact table
    event_list = clexmonte::make_complete_event_list(
        prim_event_list, prim_impact_info_list, occ_location, event_filters);

    // Construct CompleteEventCalculator
    event_calculator = std::make_shared<CompleteEventCalculator>(
        prim_event_list, prim_event_calculators, event_list.events);
  }

  /// The system
  std::shared_ptr<system_type> system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Functions for calculating event states, one for each prim event.
  /// This is supercell-specific, even though it is one per prim event,
  /// because it depends on supercell-specific clexulators
  std::vector<EventStateCalculator> prim_event_calculators;

  /// Calculator for KMC event selection
  std::shared_ptr<CompleteEventCalculator> event_calculator;
};

class KineticPotential : public BaseMontePotential {
 public:
  KineticPotential(std::shared_ptr<StateData> _state_data)
      : BaseMontePotential(_state_data),
        state(*state_data->state),
        formation_energy_clex(
            get_clex(*state_data->system, state, "formation_energy")) {}

  // --- Data used in the potential calculation: ---

  state_type const &state;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() override {
    return formation_energy_clex->per_supercell();
  }

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() override {
    return formation_energy_clex->per_unitcell();
  }

  /// \brief Calculate change in (per_supercell) potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) override {
    return formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
  }
};

class KineticCalculator : public BaseMonteCalculator {
 public:
  using BaseMonteCalculator::engine_type;

  /// Stores all the KMC event data
  std::shared_ptr<KineticEventData> event_data;

  /// Selectively allow events by unit cell
  std::vector<EventFilterGroup> event_filters;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  KineticCalculator()
      : BaseMonteCalculator("KineticCalculator",   // calculator_name
                            {},                    // required_basis_set,
                            {},                    // required_local_basis_set,
                            {"formation_energy"},  // required_clex,
                            {},                    // required_multiclex,
                            {},                    // required_local_clex,
                            {},                    // required_local_multiclex,
                            {},                    // required_dof_spaces,
                            {},                    // required_params,
                            {},                    // optional_params,
                            true,                  // time_sampling_allowed,
                            true,                  // update_atoms,
                            false,                 // save_atom_info,
                            false                  // is_multistate_method,
                            ),
        transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)) {}

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    using namespace monte_calculator;

    std::vector<state_sampling_function_type> functions =
        common_sampling_functions(
            calculation, "potential_energy",
            "Potential energy of the state (normalized per primitive cell)");

    // Specific to kmc
    functions.push_back(
        make_mean_R_squared_collective_isotropic_f(calculation));
    functions.push_back(
        make_mean_R_squared_collective_anisotropic_f(calculation));
    functions.push_back(
        make_mean_R_squared_individual_isotropic_f(calculation));
    functions.push_back(
        make_mean_R_squared_individual_anisotropic_f(calculation));
    functions.push_back(make_L_isotropic_f(calculation));
    functions.push_back(make_L_anisotropic_f(calculation));
    functions.push_back(make_D_tracer_isotropic_f(calculation));
    functions.push_back(make_D_tracer_anisotropic_f(calculation));
    functions.push_back(make_jumps_per_atom_by_type_f(calculation));
    functions.push_back(make_jumps_per_event_by_type_f(calculation));
    functions.push_back(make_jumps_per_atom_per_event_by_type_f(calculation));

    std::map<std::string, state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<json_state_sampling_function_type> functions =
        monte_calculator::common_json_sampling_functions(calculation);

    std::map<std::string, json_state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<results_analysis_function_type> functions = {
        monte_calculator::make_heat_capacity_f(calculation)};

    std::map<std::string, results_analysis_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<StateModifyingFunction> functions = {
        monte_calculator::make_match_composition_f(calculation),
        monte_calculator::make_enforce_composition_f(calculation)};

    StateModifyingFunctionMap function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct default SamplingFixtureParams
  sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results, bool write_trajectory, bool write_observations,
      bool write_status, std::optional<std::string> output_dir,
      std::optional<std::string> log_file,
      double log_frequency_in_s) const override {
    monte::SamplingParams sampling_params;
    {
      auto &s = sampling_params;
      s.sampler_names = {"clex.formation_energy",
                         "potential_energy",
                         "mol_composition",
                         "param_composition",
                         "mean_R_squared_collective_isotropic",
                         "mean_R_squared_individual_isotropic",
                         "L_isotropic",
                         "D_tracer_isotropic",
                         "mean_R_squared_collective_anisotropic",
                         "mean_R_squared_individual_anisotropic",
                         "L_anisotropic",
                         "D_tracer_anisotropic",
                         "jumps_per_atom_by_type",
                         "jumps_per_event_by_type",
                         "jumps_per_atom_per_event_by_type"};
      std::string prefix;
      prefix = "order_parameter_";
      for (auto const &pair : calculation->system()->dof_spaces) {
        s.sampler_names.push_back(prefix + pair.first);
      }
      prefix = "subspace_order_parameter_";
      for (auto const &pair : calculation->system()->dof_subspaces) {
        s.sampler_names.push_back(prefix + pair.first);
      }
      if (write_trajectory) {
        s.do_sample_trajectory = true;
      }
    }

    monte::CompletionCheckParams<statistics_type> completion_check_params;
    {
      auto &c = completion_check_params;
      c.equilibration_check_f = monte::default_equilibration_check;
      c.calc_statistics_f =
          monte::default_statistics_calculator<statistics_type>();

      completion_check_params.cutoff_params.max_count = 100;
    }

    std::vector<std::string> analysis_names = {"heat_capacity"};

    return clexmonte::make_sampling_fixture_params(
        label, calculation->sampling_functions,
        calculation->json_sampling_functions, calculation->analysis_functions,
        sampling_params, completion_check_params, analysis_names, write_results,
        write_trajectory, write_observations, write_status, output_dir,
        log_file, log_frequency_in_s);
  }

  /// \brief Validate the state's configuration
  ///
  /// Notes:
  /// - All configurations are valid (validate_state checks for consistency
  ///   with the composition conditions)
  Validator validate_configuration(state_type &state) const override {
    return Validator{};
  }

  /// \brief Validate state's conditions
  ///
  /// Notes:
  /// - requires scalar temperature
  /// - validate composition consistency
  /// - warnings if other conditions are present
  Validator validate_conditions(state_type &state) const override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in KineticCalculator::validate_conditions: system==nullptr");
    }

    // validate state.conditions
    monte::ValueMap const &conditions = state.conditions;
    Validator v;
    v.insert(validate_keys(conditions.scalar_values,
                           {"temperature"} /*required*/, {} /*optional*/,
                           "scalar", "condition", false /*throw_if_invalid*/));
    v.insert(
        validate_keys(conditions.vector_values, {} /*required*/,
                      {"param_composition", "mol_composition"} /*optional*/,
                      "vector", "condition", false /*throw_if_invalid*/));
    v.insert(validate_composition_consistency(
        state, get_composition_converter(*this->system),
        this->mol_composition_tol));
    return v;
  }

  /// \brief Validate state
  Validator validate_state(state_type &state) const override {
    Validator v;
    v.insert(this->validate_configuration(state));
    v.insert(this->validate_conditions(state));
    if (!v.valid()) {
      return v;
    }

    // check if configuration is consistent with conditions
    auto const &composition_calculator =
        get_composition_calculator(*this->system);
    auto const &composition_converter =
        get_composition_converter(*this->system);

    Eigen::VectorXd mol_composition =
        composition_calculator.mean_num_each_component(get_occupation(state));
    Eigen::VectorXd param_composition =
        composition_converter.param_composition(mol_composition);

    Eigen::VectorXd target_mol_composition =
        get_mol_composition(*this->system, state.conditions);
    Eigen::VectorXd target_param_composition =
        composition_converter.param_composition(target_mol_composition);

    if (!CASM::almost_equal(mol_composition, target_mol_composition,
                            this->mol_composition_tol)) {
      std::stringstream msg;
      msg << "***" << std::endl;
      msg << "Calculated composition is not consistent with conditions "
             "composition."
          << std::endl;
      msg << "Calculated composition:" << std::endl;
      msg << "- mol_composition: " << mol_composition.transpose() << std::endl;
      msg << "- param_composition: " << param_composition.transpose()
          << std::endl;
      msg << "Conditions:" << std::endl;
      msg << "- mol_composition: " << target_mol_composition.transpose()
          << std::endl;
      msg << "- param_composition: " << target_param_composition.transpose()
          << std::endl;
      msg << "***" << std::endl;
      v.error.insert(msg.str());
    }
    return v;
  }

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in KineticCalculator::run: system==nullptr");
    }

    // Validate state
    Validator v = this->validate_state(state);
    clexmonte::print(CASM::log(), v);
    if (!v.valid()) {
      throw std::runtime_error(
          "Error in KineticCalculator::run: Invalid initial state");
    }

    // Make state data
    this->state_data =
        std::make_shared<StateData>(this->system, &state, occ_location);

    // Make potential calculator
    this->potential = std::make_shared<KineticPotential>(this->state_data);
  }

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override {
    // Set state and potential
    // - Validates this->system is not null
    // - Validates state
    // - Throw if validation fails
    // - Constructs this->state_data
    // - Constructs this->potential
    this->set_state_and_potential(state, &occ_location);

    // if same supercell
    // -> just re-set state & avoid re-constructing event list
    if (this->transformation_matrix_to_super ==
        this->state_data->transformation_matrix_to_super) {
      for (auto &event_state_calculator :
           this->event_data->prim_event_calculators) {
        event_state_calculator.set(&state);
      }
    } else {
      this->transformation_matrix_to_super =
          this->state_data->transformation_matrix_to_super;
      this->event_data->update(state, occ_location, this->event_filters);
    }

    // Used to apply selected events: EventID -> monte::OccEvent
    auto get_event_f = [&](EventID const &selected_event_id) {
      // returns a monte::OccEvent
      return this->event_data->event_list.events.at(selected_event_id).event;
    };

    // Make selector
    // - This calculates all rates at construction
    lotto::RejectionFreeEventSelector event_selector(
        this->event_data->event_calculator,
        clexmonte::make_complete_event_id_list(
            this->state_data->n_unitcells, this->event_data->prim_event_list),
        this->event_data->event_list.impact_table,
        std::make_shared<lotto::RandomGenerator>(run_manager.engine));

    // Construct KMCData
    this->kmc_data = std::make_shared<kmc_data_type>();

    // Update atom_name_index_list
    // -- These do not change (no atoms moving to/from reservoirs) --
    auto event_system = get_event_system(*this->system);

    // Run Kinetic Monte Carlo at a single condition
    kinetic_monte_carlo_v2<EventID>(state, occ_location, *this->kmc_data,
                                    event_selector, get_event_f, run_manager,
                                    event_system);
  }

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           run_manager_type<engine_type> &run_manager) override {
    throw std::runtime_error(
        "Error: KineticCalculator does not allow multi-state runs");
  }

  // --- Parameters ---
  int verbosity_level = 10;
  double mol_composition_tol = CASM::TOL;

  /// \brief Reset the derived Monte Carlo calculator
  ///
  /// Parameters:
  ///
  ///   verbosity: str or int, default=10
  ///       If integer, the allowed range is `[0,100]`. If string, then:
  ///       - "none" is equivalent to integer value 0
  ///       - "quiet" is equivalent to integer value 5
  ///       - "standard" is equivalent to integer value 10
  ///       - "verbose" is equivalent to integer value 20
  ///       - "debug" is equivalent to integer value 100
  void _reset() override {
    ParentInputParser parser{params};

    // "verbosity": str or int, default=10
    this->verbosity_level = parse_verbosity(parser);
    CASM::log().set_verbosity(this->verbosity_level);

    // "mol_composition_tol": float, default=CASM::TOL
    this->mol_composition_tol = CASM::TOL;
    parser.optional(this->mol_composition_tol, "mol_composition_tol");

    // TODO: enumeration

    std::stringstream ss;
    ss << "Error in KineticCalculator: error reading calculation "
          "parameters.";
    std::runtime_error error_if_invalid{ss.str()};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    // Make event data  TODO: provide access to event data?
    this->event_data = std::make_shared<KineticEventData>(system);

    // TODO: Read event filters from params

    return;
  }

  /// \brief Clone the KineticCalculator
  KineticCalculator *_clone() const override {
    return new KineticCalculator(*this);
  }
};

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator() {
  return new CASM::clexmonte::kinetic_2::KineticCalculator();
}
}
