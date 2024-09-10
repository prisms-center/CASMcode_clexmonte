#ifndef CASM_clexmonte_BaseMonteCalculator
#define CASM_clexmonte_BaseMonteCalculator

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/monte_calculator/StateData.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/misc/Validator.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/ValueMap.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"
#include "casm/monte/run_management/SamplingFixture.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

/// \brief Implements a potential
class BaseMontePotential {
 public:
  BaseMontePotential(std::shared_ptr<StateData> _state_data)
      : state_data(_state_data) {}

  std::shared_ptr<StateData> state_data;

  /// \brief Calculate (per_supercell) potential value
  virtual double per_supercell() = 0;

  /// \brief Calculate (per_unitcell) potential value
  virtual double per_unitcell() = 0;

  /// \brief Calculate change in (per_supercell) semi-grand potential value due
  ///     to a series of occupation changes
  virtual double occ_delta_per_supercell(
      std::vector<Index> const &linear_site_index,
      std::vector<int> const &new_occ) = 0;
};

/// \brief Implements semi-grand canonical Monte Carlo calculations
class BaseMonteCalculator {
 public:
  typedef std::mt19937_64 engine_type;
  typedef monte::KMCData<config_type, statistics_type, engine_type>
      kmc_data_type;

  explicit BaseMonteCalculator(std::string _calculator_name,
                               std::set<std::string> _required_basis_set,
                               std::set<std::string> _required_local_basis_set,
                               std::set<std::string> _required_clex,
                               std::set<std::string> _required_multiclex,
                               std::set<std::string> _required_local_clex,
                               std::set<std::string> _required_local_multiclex,
                               std::set<std::string> _required_dof_spaces,
                               std::set<std::string> _required_params,
                               std::set<std::string> _optional_params,
                               bool _time_sampling_allowed, bool _update_atoms,
                               bool _save_atom_info,
                               bool _is_multistate_method);

  virtual ~BaseMonteCalculator();

  // --- Set at construction: ---

  /// Calculator name
  std::string calculator_name;

  /// Required system data

  std::set<std::string> required_basis_set;
  std::set<std::string> required_local_basis_set;
  std::set<std::string> required_clex;
  std::set<std::string> required_multiclex;
  std::set<std::string> required_local_clex;
  std::set<std::string> required_local_multiclex;
  std::set<std::string> required_dof_spaces;
  std::set<std::string> required_params;
  std::set<std::string> optional_params;

  /// Method allows time-based sampling?
  bool time_sampling_allowed;

  /// Method tracks species locations? (like in KMC)
  bool update_atoms;

  /// Method saves atom initial / final info? (like in KMC)
  bool save_atom_info;

  // --- Set via `reset` method: ---

  /// Calculator method parameters
  jsonParser params;

  /// System data
  std::shared_ptr<system_type> system;

  /// \brief Set parameters, check for required system data, and reset derived
  /// Monte Carlo calculator
  void reset(jsonParser const &_params, std::shared_ptr<system_type> _system) {
    this->params = _params;
    this->system = _system;
    this->_check_system();
    this->_check_params();
    this->_reset();
  }

  // --- Use after `set` and before `run` is called: ---

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  virtual std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const = 0;

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  virtual std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const = 0;

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  virtual std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const = 0;

  /// \brief Construct functions that may be used to modify states
  virtual StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const = 0;

  /// \brief Construct default SamplingFixtureParams
  virtual sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results, bool write_trajectory, bool write_observations,
      bool write_status, std::optional<std::string> output_dir,
      std::optional<std::string> log_file, double log_frequency_in_s) const = 0;

  /// \brief Validate the state's configuration
  virtual Validator validate_configuration(state_type &state) const = 0;

  /// \brief Validate the state's conditions
  virtual Validator validate_conditions(state_type &state) const = 0;

  /// \brief Validate the state
  virtual Validator validate_state(state_type &state) const = 0;

  // --- Set when `set_state_and_potential` is called: ---

  /// State data for sampling functions, for the current state
  std::shared_ptr<StateData> state_data;

  /// The current state's potential calculator, set
  ///    when the `run` method is called
  std::shared_ptr<BaseMontePotential> potential;

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  virtual void set_state_and_potential(state_type &state,
                                       monte::OccLocation *occ_location) = 0;

  // --- Set when `run` is called: ---

  /// KMC data for sampling functions, for the current state (if applicable)
  std::shared_ptr<kmc_data_type> kmc_data;

  // --- Run method: ---

  /// \brief Perform a single run, evolving current state
  virtual void run(state_type &state, monte::OccLocation &occ_location,
                   run_manager_type<engine_type> &run_manager) = 0;

  // --- Experimental, to support multi-state methods: ---

  /// Check if multi-state method
  bool is_multistate_method;

  /// Current state index
  int current_state;

  /// State data for sampling functions, to support multiple-state methods
  std::vector<std::shared_ptr<StateData>> multistate_data;

  /// State data for sampling functions, to support multiple-state methods
  std::vector<std::shared_ptr<BaseMontePotential>> multistate_potential;

  /// \brief Perform a single run, evolving one or more states
  virtual void run(int current_state, std::vector<state_type> &states,
                   std::vector<monte::OccLocation> &occ_locations,
                   run_manager_type<engine_type> &run_manager) = 0;

  /// \brief Clone the BaseMonteCalculator
  std::unique_ptr<BaseMonteCalculator> clone() const;

 private:
  /// \brief Standardized check for whether system has required data
  void _check_system() const;

  /// \brief Validate input JSON (for key existence only)
  void _check_params() const;

  /// \brief Reset the derived Monte Carlo calculator
  virtual void _reset() = 0;

  /// \brief Clone the BaseMonteCalculator
  virtual BaseMonteCalculator *_clone() const = 0;
};

/// \brief Validate input parameters (for key existence only)
template <typename T>
Validator validate_keys(std::map<std::string, T> const &map,
                        std::set<std::string> required,
                        std::set<std::string> optional, std::string which_type,
                        std::string kind, bool throw_if_invalid = true);

/// \brief Print validation results
void print(Log &log, Validator const &validator);

// ~~~ Implementations ---

/// \brief Validate input parameters (for key existence only)
///
/// Notes:
/// - Input parameters that start with "_" are ignored
/// - A warning is printed for input parameters that are not in the
///   required or optional sets
///
/// \tparam T Input parameter type
/// \param map Input parameter container as key:value pairs
/// \param required Required input parameters
/// \param optional Optional input parameters
/// \param which_type One of "bool", "float", "vector", "matrix", "str"
/// \param kind One of "parameter", "condition", etc.
template <typename T>
Validator validate_keys(std::map<std::string, T> const &map,
                        std::set<std::string> required,
                        std::set<std::string> optional, std::string which_type,
                        std::string kind, bool throw_if_invalid) {
  Validator validator;
  for (auto key : required) {
    if (!map.count(key)) {
      std::stringstream msg;
      msg << "Error: Missing required " << which_type << " " << kind << " '"
          << key << "'.";
      if (throw_if_invalid) {
        throw std::runtime_error(msg.str());
      } else {
        validator.error.insert(msg.str());
      }
    }
  }

  auto &log = CASM::log();
  for (auto const &pair : map) {
    std::string key = pair.first;
    if (key.empty()) {
      std::stringstream msg;
      msg << "Error: Empty " << which_type << " " << kind << " value.";
      if (throw_if_invalid) {
        throw std::runtime_error(msg.str());
      } else {
        validator.error.insert(msg.str());
      }
    }
    if (key[0] == '_') {
      continue;
    }
    if (!required.count(key) && !optional.count(key)) {
      std::stringstream msg;
      msg << "Warning: Unknown " << which_type << " " << kind << " '" << key
          << "'.";
      if (throw_if_invalid) {
        log.indent() << msg.str() << std::endl;
      } else {
        validator.warning.insert(msg.str());
      }
    }
  }
  return validator;
}

inline void print(Log &log, Validator const &validator) {
  if (!validator.valid()) {
    log.custom("Errors");
    for (auto const &msg : validator.error) {
      log.indent() << "- " << msg << std::endl;
    }
    log << std::endl;
  }
  if (validator.warning.size()) {
    log.custom("Warnings");
    for (auto const &msg : validator.warning) {
      log.indent() << "- " << msg << std::endl;
    }
    log << std::endl;
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
