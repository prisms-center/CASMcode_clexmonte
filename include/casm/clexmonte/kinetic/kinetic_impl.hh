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
#include "casm/monte/results/Results.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

// debug
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/events/io/json/EventState_json_io.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Implements kinetic Monte Carlo calculations
template <typename EngineType>
Kinetic<EngineType>::Kinetic(std::shared_ptr<system_type> _system,
                             std::shared_ptr<EngineType> _random_number_engine)
    : system(_system),
      random_number_generator(_random_number_engine),
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
                              run_manager_type &run_manager) {
  this->state = &state;
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);
  Index n_unitcells = this->transformation_matrix_to_super.determinant();

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
    this->event_data->update(state, this->conditions, occ_location);
  }

  // Enforce composition -- occ_location is maintained up-to-date
  monte::Conversions const &convert = get_index_conversions(*system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*system, state);
  std::vector<monte::OccSwap> grand_canonical_swaps =
      make_grand_canonical_swaps(convert, occ_candidate_list);
  clexmonte::enforce_composition(
      get_occupation(state),
      state.conditions.vector_values.at("mol_composition"),
      get_composition_calculator(*system), grand_canonical_swaps, occ_location,
      this->random_number_generator);

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
      this->event_data->event_list.impact_table);

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
      make_kmc_potential_energy_f(calculation),
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

  std::map<std::string, state_sampling_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
};

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
std::map<std::string, state_modifying_function_type>
Kinetic<EngineType>::standard_modifying_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  std::vector<state_modifying_function_type> functions = {
      make_set_mol_composition_f(calculation)};

  std::map<std::string, state_modifying_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
