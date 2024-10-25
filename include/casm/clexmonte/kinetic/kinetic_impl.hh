#ifndef CASM_clexmonte_kinetic_impl
#define CASM_clexmonte_kinetic_impl

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/kinetic/kinetic.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/kinetic_sampling_functions.hh"
#include "casm/clexmonte/state/modifying_functions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/run_management/Results.hh"
#include "casm/monte/run_management/State.hh"

// debug
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Implements kinetic Monte Carlo calculations
template <typename EngineType>
Kinetic<EngineType>::Kinetic(std::shared_ptr<system_type> _system,
                             std::vector<EventFilterGroup> _event_filters)
    : system(_system),
      event_filters(_event_filters),
      event_data(std::make_shared<KineticEventData>(system)),
      state(nullptr),
      transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)),
      occ_location(nullptr) {
  if (!is_clex_data(*this->system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing Kinetic: no 'formation_energy' clex.");
  }
}

/// \brief Perform a single run, evolving current state
template <typename EngineType>
void Kinetic<EngineType>::run(state_type &state,
                              monte::OccLocation &occ_location,
                              run_manager_type<EngineType> &run_manager) {
  this->state = &state;
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);
  Index n_unitcells = this->transformation_matrix_to_super.determinant();

  // Make potential calculator - for sampling function only
  this->potential =
      std::make_shared<canonical::CanonicalPotential>(this->system);
  this->potential->set(this->state, this->conditions);
  this->formation_energy = this->potential->formation_energy();

  // if same supercell
  // -> just re-set state & conditions & avoid re-constructing event list
  if (this->transformation_matrix_to_super ==
          get_transformation_matrix_to_super(state) &&
      this->conditions != nullptr) {
    for (auto &event_state_calculator :
         this->event_data->prim_event_calculators) {
      event_state_calculator.set(this->state, this->conditions);
    }
  } else {
    this->transformation_matrix_to_super =
        get_transformation_matrix_to_super(state);
    n_unitcells = this->transformation_matrix_to_super.determinant();
    this->event_data->update(state, this->conditions, occ_location,
                             this->event_filters);
  }

  // Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator(
      run_manager.engine);

  // Enforce composition -- occ_location is maintained up-to-date
  std::vector<monte::OccSwap> const &semigrand_canonical_swaps =
      get_semigrand_canonical_swaps(*this->system);
  clexmonte::enforce_composition(
      get_occupation(state),
      state.conditions.vector_values.at("mol_composition"),
      get_composition_calculator(*system), semigrand_canonical_swaps,
      occ_location, random_number_generator);

  // Used to apply selected events: EventID -> monte::OccEvent
  auto get_event_f = [&](EventID const &selected_event_id) {
    // returns a monte::OccEvent
    return this->event_data->event_list.events.at(selected_event_id).event;
  };

  // Make selector
  lotto::RejectionFreeEventSelector event_selector(
      this->event_data->event_calculator,
      clexmonte::make_complete_event_id_list(n_unitcells,
                                             this->event_data->prim_event_list),
      this->event_data->event_list.impact_table,
      std::make_shared<lotto::RandomGenerator>(run_manager.engine));

  // Update atom_name_index_list -- These do not change --
  // TODO: KMC with atoms that move to/from resevoir will need to update this
  auto event_system = get_event_system(*this->system);
  this->kmc_data.atom_name_index_list =
      make_atom_name_index_list(occ_location, *event_system);

  monte::kinetic_monte_carlo<EventID>(state, occ_location, this->kmc_data,
                                      event_selector, get_event_f, run_manager);
}

/// \brief Construct functions that may be used to sample various quantities
///     of the Monte Carlo calculation as it runs
template <typename EngineType>
std::map<std::string, state_sampling_function_type>
Kinetic<EngineType>::standard_sampling_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  std::vector<state_sampling_function_type> functions = {
      make_temperature_f(calculation),
      make_mol_composition_f(calculation),
      make_param_composition_f(calculation),
      make_formation_energy_corr_f(calculation),
      make_formation_energy_f(calculation),
      make_potential_energy_f(calculation),
      make_mean_R_squared_collective_isotropic_f(calculation),
      make_mean_R_squared_collective_anisotropic_f(calculation),
      make_mean_R_squared_individual_isotropic_f(calculation),
      make_mean_R_squared_individual_anisotropic_f(calculation),
      make_L_isotropic_f(calculation),
      make_L_anisotropic_f(calculation),
      make_D_tracer_isotropic_f(calculation),
      make_D_tracer_anisotropic_f(calculation),
      make_jumps_per_atom_by_type_f(calculation),
      make_jumps_per_event_by_type_f(calculation),
      make_jumps_per_atom_per_event_by_type_f(calculation)};

  make_order_parameter_f(functions, calculation);
  make_subspace_order_parameter_f(functions, calculation);

  std::map<std::string, state_sampling_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param calculation Shared pointer to Canonical calculation, which
///     can be used by sampling functions to access system and calculation data
///     such as the prim, the cluster expansion, and the composition axes.
///
template <typename EngineType>
std::map<std::string, json_state_sampling_function_type>
Kinetic<EngineType>::standard_json_sampling_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  std::vector<json_state_sampling_function_type> functions = {
      make_config_f(calculation)};

  std::map<std::string, json_state_sampling_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
template <typename EngineType>
std::map<std::string, results_analysis_function_type>
Kinetic<EngineType>::standard_analysis_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  std::vector<results_analysis_function_type> functions = {
      make_heat_capacity_f(calculation), make_mol_susc_f(calculation),
      make_param_susc_f(calculation), make_mol_thermochem_susc_f(calculation),
      make_param_thermochem_susc_f(calculation)};

  std::map<std::string, results_analysis_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to modify states
template <typename EngineType>
StateModifyingFunctionMap Kinetic<EngineType>::standard_modifying_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  std::vector<StateModifyingFunction> functions = {
      make_set_mol_composition_f(calculation)};

  StateModifyingFunctionMap function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
