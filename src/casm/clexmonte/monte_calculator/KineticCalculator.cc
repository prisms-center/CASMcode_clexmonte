#include "casm/clexmonte/monte_calculator/KineticCalculator.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/methods/kinetic_monte_carlo.hh"
#include "casm/clexmonte/monte_calculator/MonteEventData.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/kinetic_events.hh"
#include "casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/modifying_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/monte_calculator/selected_event_data_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"
#include "casm/monte/sampling/io/json/SelectedEventData_json_io.hh"

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
                          {"verbosity", "mol_composition_tol",
                           "event_data_type", "event_selector_type",
                           "impact_table_type", "allow_events_with_no_barrier",
                           "selected_event_data"},  // optional_params,
                          true,                     // time_sampling_allowed,
                          true,                     // update_atoms,
                          false,                    // save_atom_info,
                          false                     // is_multistate_method,
      ) {
  // this could go into base constructor
  this->selected_event = std::make_shared<SelectedEvent>();
  this->selected_event_data_functions =
      std::make_shared<monte::SelectedEventDataFunctions>();
  this->selected_event_data = std::make_shared<monte::SelectedEventData>();
}

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
  functions.push_back(make_selected_event_count_by_type_f(calculation));
  functions.push_back(make_selected_event_fraction_by_type_f(calculation));
  functions.push_back(
      make_selected_event_count_by_equivalent_index_f(calculation));
  functions.push_back(
      make_selected_event_fraction_by_equivalent_index_f(calculation));
  functions.push_back(
      make_selected_event_count_by_equivalent_index_and_direction_f(
          calculation));
  functions.push_back(
      make_selected_event_fraction_by_equivalent_index_and_direction_f(
          calculation));

  for (auto f : make_selected_event_count_by_equivalent_index_per_event_type_f(
           calculation)) {
    functions.push_back(f);
  }
  for (auto f :
       make_selected_event_fraction_by_equivalent_index_per_event_type_f(
           calculation)) {
    functions.push_back(f);
  }

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

/// \brief Construct functions that may be used to collect selected event data
std::optional<monte::SelectedEventDataFunctions>
KineticCalculator::standard_selected_event_data_functions(
    std::shared_ptr<MonteCalculator> const &calculation) const {
  using namespace monte_calculator;
  monte::SelectedEventDataFunctions functions;

  // Event type data:
  functions.insert(make_selected_event_by_type_f(calculation));
  functions.insert(make_selected_event_by_equivalent_index_f(calculation));
  functions.insert(
      make_selected_event_by_equivalent_index_and_direction_f(calculation));
  for (auto f :
       make_selected_event_by_equivalent_index_per_event_type_f(calculation)) {
    functions.insert(f);
  }
  for (auto f : make_local_orbit_composition_f(calculation)) {
    functions.insert(f);
  }

  // Event state data:
  functions.insert(make_dE_activated_by_type_f(calculation));
  functions.insert(make_dE_activated_by_equivalent_index_f(calculation));

  return functions;
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
    s.do_sample_time = true;
    std::string prefix;
    prefix = "order_parameter.";
    for (auto const &pair : calculation->system()->dof_spaces) {
      s.sampler_names.push_back(prefix + pair.first);
    }
    prefix = "order_parameter.";
    std::string suffix = ".subspace_magnitudes";
    for (auto const &pair : calculation->system()->dof_subspaces) {
      s.sampler_names.push_back(prefix + pair.first + suffix);
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
///
/// \param state State to set
/// \param occ_location Pointer to OccLocation to use, or may be nullptr
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
/// Notes:
/// - Validates this->state_data is not null
/// - Validates this->state_data->occ_location is not null
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

  // Currently, event_filters are only set at _reset() by reading from params
  std::optional<std::vector<EventFilterGroup>> event_filters = std::nullopt;
  if (this->event_data_type == kinetic_event_data_type::high_memory) {
    this->_complete_event_data().update(state, occ_location, event_filters,
                                        engine);
  } else if (this->event_data_type == kinetic_event_data_type::default_memory) {
    this->_allowed_event_data().update(state, occ_location, engine);
  } else if (this->event_data_type == kinetic_event_data_type::low_memory) {
    this->_allowed_event_data().update(state, occ_location, engine);

  } else {
    throw std::runtime_error(
        "Error in KineticCalculator::set_event_data: event_data_type is "
        "invalid");
  }
}

/// \brief Perform a single run, evolving current state
void KineticCalculator::run(state_type &state, monte::OccLocation &occ_location,
                            run_manager_type<engine_type> &run_manager) {
  if (run_manager.sampling_fixtures.size() == 0) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: "
        "run_manager.sampling_fixtures.size()==0");
  }

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
  std::cout << "Setting event data ... " << std::endl;
  this->set_event_data(run_manager.engine);
  std::cout << "Setting event data ... DONE" << std::endl << std::endl;

  // Construct EventDataSummary
  std::cout << "Generating event data summary ... " << std::flush;
  MonteEventData monte_event_data(this->event_data, nullptr);
  double energy_bin_width = 0.1;
  double freq_bin_width = 0.1;
  double rate_bin_width = 0.1;

  EventDataSummary event_data_summary(this->state_data, monte_event_data,
                                      energy_bin_width, freq_bin_width,
                                      rate_bin_width);
  std::cout << "Generating event data summary ... DONE" << std::endl
            << std::endl;
  print(std::cout, event_data_summary);

  // Construct KMCData
  this->kmc_data = std::make_shared<kmc_data_type>();

  // Update atom_name_index_list
  // -- These do not change (no atoms moving to/from reservoirs) --
  auto event_system = get_event_system(*this->system);

  // Optional: Manages constructing histogram data structures
  // and collecting selected event data

  std::optional<monte::SelectedEventDataCollector> collector;
  if (this->selected_event_data_params) {
    if (this->selected_event_data_functions == nullptr) {
      throw std::runtime_error(
          "Error in KineticCalculator::run: "
          "this->selected_event_data_functions==nullptr");
    }
    collector = monte::SelectedEventDataCollector(
        *selected_event_data_functions, *selected_event_data_params,
        selected_event_data);
  }

  // Check this->selected_event is not null
  if (this->selected_event == nullptr) {
    throw std::runtime_error(
        "Error in KineticCalculator::run: this->selected_event==nullptr");
  }

  // Function to set selected event
  bool requires_event_state =
      collector.has_value() && collector->requires_event_state;
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->event_data->select_event(selected_event, requires_event_state);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<EventID>(state, occ_location, *this->kmc_data,
                                  *this->selected_event, set_selected_event_f,
                                  collector, run_manager, event_system);

  // Warn if events with no barrier were encountered
  check_n_not_normal(this->event_data->n_not_normal());
}

/// \brief Print a warning to std::cerr if events with no barrier were
///     encountered
void KineticCalculator::check_n_not_normal(
    std::map<std::string, Index> const &n_not_normal) const {
  if (n_not_normal.empty()) {
    return;
  }

  std::cerr << "## WARNING: EVENTS WITH NO BARRIER ##################\n"
               "#                                                   #\n"
               "# Events with no barrier are treated as having a    #\n"
               "# rate equal to the attempt frequency.              #\n"
               "#                                                   #\n"
               "# Number encountered by type:                       #\n";
  for (auto const &pair : n_not_normal) {
    std::cerr << "  - " << pair.first << ": " << pair.second << "\n";
  }
  std::cerr << "#                                                   #\n"
               "#####################################################\n"
            << std::endl;
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
  // -- Parsing ----------------------------

  // `params` is BaseMonteCalculator::params
  // - originally "calculation_options" from the run params input file
  // - `params` argument of MonteCalculator Python constructor

  ParentInputParser parser{params};

  // "verbosity": str or int, default=10
  this->verbosity_level = parse_verbosity(parser);
  CASM::log().set_verbosity(this->verbosity_level);

  // "mol_composition_tol": float, default=CASM::TOL
  this->mol_composition_tol = CASM::TOL;
  parser.optional(this->mol_composition_tol, "mol_composition_tol");

  // TODO: enumeration

  // TODO: Read event_filters from params
  std::optional<std::vector<EventFilterGroup>> event_filters = std::nullopt;

  // Read selected event data params
  this->selected_event_data_params.reset();
  if (parser.self.contains("selected_event_data")) {
    auto selected_event_data_subparser =
        parser.subparse<monte::SelectedEventDataParams>("selected_event_data");
    if (selected_event_data_subparser->valid()) {
      this->selected_event_data_params =
          std::move(selected_event_data_subparser->value);
    }
  }

  // Read "event_data_type"
  // - "high_memory": complete event list,
  // - "default" (default): allowed event list w/ vector index, or
  // - "low_memory": allowed event list w/ map index
  std::string event_data_type_str = "default";
  parser.optional(event_data_type_str, "event_data_type");
  if (event_data_type_str == "high_memory") {
    this->event_data_type = kinetic_event_data_type::high_memory;
  } else if (event_data_type_str == "default") {
    this->event_data_type = kinetic_event_data_type::default_memory;
  } else if (event_data_type_str == "low_memory") {
    this->event_data_type = kinetic_event_data_type::low_memory;

    if (event_filters.has_value()) {
      parser.insert_error(
          "event_data_type",
          "event_filters are not supported by event_data_type 'low_memory'");
    }
  } else {
    parser.insert_error("event_data_type",
                        "Invalid event_data_type: " + event_data_type_str);
  }

  // Read "event_selector_type"
  // - "vector_sum_tree" (default): binary sum tree (log complexity) built
  //   using std::vector
  // - "sum_tree": binary sum tree (log complexity) built using std::map and
  //   linked lists
  // - "direct_sum": direct sum of rate vector (linear complexity)
  std::string event_selector_type_str = "vector_sum_tree";
  parser.optional(event_selector_type_str, "event_selector_type");
  if (event_selector_type_str == "vector_sum_tree") {
    this->event_selector_type = kinetic_event_selector_type::vector_sum_tree;
  } else if (event_selector_type_str == "sum_tree") {
    this->event_selector_type = kinetic_event_selector_type::sum_tree;
  } else if (event_selector_type_str == "direct_sum") {
    this->event_selector_type = kinetic_event_selector_type::direct_sum;
  } else {
    parser.insert_error("event_selector_type", "Invalid event_selector_type: " +
                                                   event_selector_type_str);
  }

  /// Read "impact_table_type" (optional, only takes affect for "low_memory")
  this->use_neighborlist_impact_table = true;
  std::string impact_table_type_str = "neighborlist";
  parser.optional(impact_table_type_str, "impact_table_type");
  if (impact_table_type_str == "neighborlist") {
    this->use_neighborlist_impact_table = true;
  } else if (impact_table_type_str == "relative") {
    this->use_neighborlist_impact_table = false;
  } else {
    parser.insert_error("impact_table_type",
                        "Invalid impact_table_type: " + impact_table_type_str);
  }

  // Read "allow_events_with_no_barrier"
  // - true or false (default)
  this->allow_events_with_no_barrier = false;
  parser.optional(allow_events_with_no_barrier, "allow_events_with_no_barrier");

  std::stringstream ss;
  ss << "Error in KineticCalculator: error reading calculation "
        "parameters.";
  std::runtime_error error_if_invalid{ss.str()};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  // -- After parsing ----------------------------

  // Make event data
  if (this->event_data_type == kinetic_event_data_type::high_memory) {
    this->event_data = std::make_shared<kinetic_2::CompleteKineticEventData>(
        system, event_filters, this->allow_events_with_no_barrier);
  } else if (this->event_data_type == kinetic_event_data_type::default_memory) {
    this->event_data = std::make_shared<kinetic_2::AllowedKineticEventData>(
        system, this->allow_events_with_no_barrier, false /* use_map_index */,
        this->use_neighborlist_impact_table, this->event_selector_type);
  } else {
    this->event_data = std::make_shared<kinetic_2::AllowedKineticEventData>(
        system, this->allow_events_with_no_barrier, true /* use_map_index */,
        this->use_neighborlist_impact_table, this->event_selector_type);
  }

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
