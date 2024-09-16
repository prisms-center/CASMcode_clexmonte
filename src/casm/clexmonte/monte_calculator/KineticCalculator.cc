#include "casm/clexmonte/monte_calculator/KineticCalculator.hh"

#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/kinetic_events.hh"
#include "casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/modifying_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

KineticPotential::KineticPotential(std::shared_ptr<StateData> _state_data)
    : BaseMontePotential(_state_data),
      state(*state_data->state),
      formation_energy_clex(
          get_clex(*state_data->system, state, "formation_energy")) {}

/// \brief Calculate (per_supercell) potential value
double KineticPotential::per_supercell() {
  return formation_energy_clex->per_supercell();
}

/// \brief Calculate (per_unitcell) potential value
double KineticPotential::per_unitcell() {
  return formation_energy_clex->per_unitcell();
}

/// \brief Calculate change in (per_supercell) potential value due
///     to a series of occupation changes
double KineticPotential::occ_delta_per_supercell(
    std::vector<Index> const &linear_site_index,
    std::vector<int> const &new_occ) {
  return formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
}

KineticCalculator::KineticCalculator()
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
KineticCalculator::standard_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  using namespace monte_calculator;

  std::vector<state_sampling_function_type> functions =
      common_sampling_functions(
          calculation, "potential_energy",
          "Potential energy of the state (normalized per primitive cell)");

  // Specific to kmc
  functions.push_back(make_mean_R_squared_collective_isotropic_f(calculation));
  functions.push_back(
      make_mean_R_squared_collective_anisotropic_f(calculation));
  functions.push_back(make_mean_R_squared_individual_isotropic_f(calculation));
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
KineticCalculator::standard_json_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
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
KineticCalculator::standard_analysis_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  std::vector<results_analysis_function_type> functions = {
      monte_calculator::make_heat_capacity_f(calculation)};

  std::map<std::string, results_analysis_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to modify states
StateModifyingFunctionMap KineticCalculator::standard_modifying_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
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
sampling_fixture_params_type
KineticCalculator::make_default_sampling_fixture_params(
    std::shared_ptr<MonteCalculator> const &calculation, std::string label,
    bool write_results, bool write_trajectory, bool write_observations,
    bool write_status, std::optional<std::string> output_dir,
    std::optional<std::string> log_file, double log_frequency_in_s) const {
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
      write_trajectory, write_observations, write_status, output_dir, log_file,
      log_frequency_in_s);
}

/// \brief Validate the state's configuration
///
/// Notes:
/// - All configurations are valid (validate_state checks for consistency
///   with the composition conditions)
Validator KineticCalculator::validate_configuration(state_type &state) const {
  return Validator{};
}

/// \brief Validate state's conditions
///
/// Notes:
/// - requires scalar temperature
/// - validate composition consistency
/// - warnings if other conditions are present
Validator KineticCalculator::validate_conditions(state_type &state) const {
  // Validate system
  if (this->system == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::validate_conditions: system==nullptr");
  }

  // validate state.conditions
  monte::ValueMap const &conditions = state.conditions;
  Validator v;
  v.insert(validate_keys(conditions.scalar_values, {"temperature"} /*required*/,
                         {} /*optional*/, "scalar", "condition",
                         false /*throw_if_invalid*/));
  v.insert(validate_keys(conditions.vector_values, {} /*required*/,
                         {"param_composition", "mol_composition"} /*optional*/,
                         "vector", "condition", false /*throw_if_invalid*/));
  v.insert(validate_composition_consistency(
      state, get_composition_converter(*this->system),
      this->mol_composition_tol));
  return v;
}

/// \brief Validate state
Validator KineticCalculator::validate_state(state_type &state) const {
  Validator v;
  v.insert(this->validate_configuration(state));
  v.insert(this->validate_conditions(state));
  if (!v.valid()) {
    return v;
  }

  // check if configuration is consistent with conditions
  auto const &composition_calculator =
      get_composition_calculator(*this->system);
  auto const &composition_converter = get_composition_converter(*this->system);

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
void KineticCalculator::set_state_and_potential(
    state_type &state, monte::OccLocation *occ_location) {
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

/// \brief Set event data (includes calculating all rates), using current
/// state data
///
/// \param engine The random number generator engine used to select events
///     and timesteps. If nullptr, a new engine is constructed,
///     seeded by std::random_device
void KineticCalculator::set_event_data(std::shared_ptr<engine_type> engine) {
  if (this->state_data == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: "
        "this->state_data==nullptr");
  }
  if (this->state_data->occ_location == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: "
        "this->state_data->occ_location==nullptr");
  }

  state_type const &state = *this->state_data->state;
  monte::OccLocation const &occ_location = *this->state_data->occ_location;

  // if same supercell
  // -> just re-set state & avoid re-constructing event list
  if (this->transformation_matrix_to_super ==
      this->state_data->transformation_matrix_to_super) {
    for (auto &event_state_calculator : _event_data().prim_event_calculators) {
      event_state_calculator.set(&state);
    }
  } else {
    this->transformation_matrix_to_super =
        this->state_data->transformation_matrix_to_super;
    _event_data().update(state, occ_location, this->event_filters);
  }

  // Make event selector
  // - This calculates all rates at construction
  _event_data().event_selector =
      std::make_shared<KineticEventData::event_selector_type>(
          _event_data().event_calculator,
          clexmonte::make_complete_event_id_list(this->state_data->n_unitcells,
                                                 _event_data().prim_event_list),
          _event_data().event_list.impact_table,
          std::make_shared<lotto::RandomGenerator>(engine));
}

/// \brief Perform a single run, evolving current state
void KineticCalculator::run(state_type &state, monte::OccLocation &occ_location,
                            run_manager_type<engine_type> &run_manager) {
  // Set state and potential
  // - Validates this->system is not null
  // - Validates state
  // - Throw if validation fails
  // - Constructs this->state_data
  // - Constructs this->potential
  this->set_state_and_potential(state, &occ_location);

  // Set event data
  // - Set or re-set state in event_state_calculators
  // - Update event_data if supercell has changed
  // - Throw if this->state_data is null
  // - Constructs this->event_data->event_selector
  // - Calculates all rates
  this->set_event_data(run_manager.engine);

  // Used to apply selected events: EventID -> monte::OccEvent
  auto get_event_f = [&](EventID const &selected_event_id) {
    // returns a monte::OccEvent
    return _event_data().event_list.events.at(selected_event_id).event;
  };

  // Construct KMCData
  this->kmc_data = std::make_shared<kmc_data_type>();

  // Update atom_name_index_list
  // -- These do not change (no atoms moving to/from reservoirs) --
  auto event_system = get_event_system(*this->system);

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<EventID>(state, occ_location, *this->kmc_data,
                                  *_event_data().event_selector, get_event_f,
                                  run_manager, event_system);
}

/// \brief Perform a single run, evolving one or more states
void KineticCalculator::run(int current_state, std::vector<state_type> &states,
                            std::vector<monte::OccLocation> &occ_locations,
                            run_manager_type<engine_type> &run_manager) {
  throw std::runtime_error(
      "Error: KineticCalculator does not allow multi-state runs");
}

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
void KineticCalculator::_reset() {
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
KineticCalculator *KineticCalculator::_clone() const {
  return new KineticCalculator(*this);
}

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator() {
  return new CASM::clexmonte::kinetic_2::KineticCalculator();
}
}
