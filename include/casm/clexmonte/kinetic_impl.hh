#ifndef CASM_clexmonte_kinetic_impl
#define CASM_clexmonte_kinetic_impl

#include <random>

#include "casm/clexmonte/canonical.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/kinetic.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/kinetic_sampling_functions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Implements kinetic Monte Carlo calculations
template <typename EngineType>
Kinetic<EngineType>::Kinetic(std::shared_ptr<system_type> _system,
                             std::shared_ptr<EngineType> _random_number_engine)
    : system(_system),
      random_number_generator(_random_number_engine),
      state(nullptr),
      transformation_matrix_to_supercell(Eigen::Matrix3l::Zero(3, 3)),
      occ_location(nullptr) {
  if (!is_clex_data(*this->system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing Kinetic: no 'formation_energy' clex.");
  }

  this->prim_event_list = clexmonte::make_prim_event_list(*this->system);

  this->prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *this->system, this->prim_event_list, {"formation_energy"});
}

/// \brief Perform a single run, evolving current state
template <typename EngineType>
void Kinetic<EngineType>::run(state_type &state,
                              monte::OccLocation &occ_location,
                              monte::RunManager<config_type> &run_manager) {
  // if same supercell
  // -> just re-set state & conditions & avoid re-constructing event list
  if (this->transformation_matrix_to_supercell ==
          get_transformation_matrix_to_supercell(state) &&
      this->conditions != nullptr) {
    this->state = &state;
    this->occ_location = &occ_location;
    *this->conditions = *make_conditions(*this->system, state);

    for (auto &event_state_calculator : this->prim_event_calculators) {
      event_state_calculator.set(this->state, this->conditions);
    }
  } else {
    this->state = &state;
    this->transformation_matrix_to_supercell =
        get_transformation_matrix_to_supercell(state);
    this->occ_location = &occ_location;
    this->conditions = make_conditions(*this->system, state);

    // These are constructed/re-constructed so cluster expansions point
    // at the current state
    this->prim_event_calculators = clexmonte::make_prim_event_calculators(
        this->system, state, this->prim_event_list, this->conditions);

    // TODO: rejection-kmc option does not require impact table
    this->event_list = clexmonte::make_complete_event_list(
        this->prim_event_list, this->prim_impact_info_list, occ_location);

    // Construct CompleteEventCalculator
    this->event_calculator =
        std::make_shared<clexmonte::CompleteEventCalculator>(
            this->prim_event_list, this->prim_event_calculators,
            this->event_list.events);
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
    return this->event_list.events.at(selected_event_id).event;
  };

  // Make selector
  Eigen::Matrix3l T = get_transformation_matrix_to_supercell(state);
  lotto::RejectionFreeEventSelector event_selector(
      this->event_calculator,
      clexmonte::make_complete_event_id_list(T.determinant(),
                                             this->prim_event_list),
      this->event_list.impact_table);

  monte::kinetic_monte_carlo<EventID>(state, occ_location, event_selector,
                                      get_event_f, run_manager);
}

/// \brief Perform a series of runs, according to a state generator
template <typename EngineType>
void Kinetic<EngineType>::run_series(
    state_generator_type &state_generator,
    std::vector<monte::SamplingFixtureParams<config_type>> const
        &sampling_fixture_params) {
  clexmonte::run_series(*this, state_generator, sampling_fixture_params);
}

/// \brief Construct functions that may be used to sample various quantities
///     of the Monte Carlo calculation as it runs
template <typename EngineType>
monte::StateSamplingFunctionMap<Configuration>
Kinetic<EngineType>::standard_sampling_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  auto const &system = calculation->system;
  std::vector<monte::StateSamplingFunction<Configuration>> functions = {
      make_temperature_f(system),
      make_mol_composition_f(system),
      make_param_composition_f(system),
      make_formation_energy_corr_f(system),
      make_formation_energy_f(system),
      make_kmc_potential_energy_f(calculation),
      make_R_squared_center_f(calculation),
      make_R_squared_tracer_f(calculation)};

  monte::StateSamplingFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
};

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
template <typename EngineType>
monte::ResultsAnalysisFunctionMap<Configuration>
Kinetic<EngineType>::standard_analysis_functions(
    std::shared_ptr<Kinetic<EngineType>> const &calculation) {
  auto const &system = calculation->system;
  std::vector<monte::ResultsAnalysisFunction<Configuration>> functions = {
      make_heat_capacity_f(), make_mol_susc_f(system),
      make_param_susc_f(system), make_mol_thermochem_susc_f(system),
      make_param_thermochem_susc_f(system)};

  monte::ResultsAnalysisFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
