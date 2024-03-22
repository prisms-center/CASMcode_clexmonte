#ifndef CASM_clexmonte_semigrand_canonical_impl
#define CASM_clexmonte_semigrand_canonical_impl

#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/methods/occupation_metropolis.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/semigrand_canonical/calculator.hh"
#include "casm/clexmonte/semigrand_canonical/event_generator.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/modifying_functions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/run_management/Results.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

template <typename EngineType>
SemiGrandCanonical<EngineType>::SemiGrandCanonical(
    std::shared_ptr<system_type> _system)
    : system(_system),
      state(nullptr),
      transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)),
      occ_location(nullptr) {
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing SemiGrandCanonical: no 'formation_energy' clex.");
  }
}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void SemiGrandCanonical<EngineType>::run(
    state_type &state, monte::OccLocation &occ_location,
    run_manager_type<EngineType> &run_manager) {
  // Store state data, which makes it available to samplers
  this->state = &state;
  this->transformation_matrix_to_super =
      get_transformation_matrix_to_super(state);
  this->occ_location = &occ_location;
  this->conditions = std::make_shared<SemiGrandCanonicalConditions>(
      get_composition_converter(*this->system));
  this->conditions->set_all(state.conditions, false);

  // Make potential calculator
  this->potential = std::make_shared<SemiGrandCanonicalPotential>(this->system);
  this->potential->set(this->state, this->conditions);
  this->formation_energy = this->potential->formation_energy();

  auto potential_occ_delta_per_supercell_f = [=](monte::OccEvent const &event) {
    return this->potential->occ_delta_per_supercell(event.linear_site_index,
                                                    event.new_occ);
  };

  // Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator(
      run_manager.engine);

  // Make event generator
  auto event_generator =
      std::make_shared<SemiGrandCanonicalEventGenerator<EngineType>>(
          get_semigrand_canonical_swaps(*this->system),
          get_semigrand_canonical_multiswaps(*this->system));
  event_generator->set(&state, &occ_location);

  auto propose_event_f =
      [=](monte::RandomNumberGenerator<engine_type> &random_number_generator)
      -> monte::OccEvent const & {
    return event_generator->propose(random_number_generator);
  };

  auto apply_event_f = [=](monte::OccEvent const &occ_event) -> void {
    return event_generator->apply(occ_event);
  };

  // Run Monte Carlo at a single condition
  clexmonte::occupation_metropolis_v2(
      state, occ_location, this->conditions->temperature,
      potential_occ_delta_per_supercell_f, propose_event_f, apply_event_f,
      run_manager);
}

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param calculation Shared pointer to SemiGrandCanonical calculation, which
///     can be used by sampling functions to access system data, such as the
///     prim, the cluster expansion, and the composition axes, and calculation
///     data, such as the potential.
///
template <typename EngineType>
std::map<std::string, state_sampling_function_type>
SemiGrandCanonical<EngineType>::standard_sampling_functions(
    std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation) {
  std::vector<state_sampling_function_type> functions = {
      make_temperature_f(calculation),
      make_mol_composition_f(calculation),
      make_param_composition_f(calculation),
      make_param_chem_pot_f(calculation),
      make_formation_energy_corr_f(calculation),
      make_formation_energy_f(calculation),
      make_potential_energy_f(calculation)};

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
SemiGrandCanonical<EngineType>::standard_json_sampling_functions(
    std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation) {
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
SemiGrandCanonical<EngineType>::standard_analysis_functions(
    std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation) {
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
StateModifyingFunctionMap
SemiGrandCanonical<EngineType>::standard_modifying_functions(
    std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation) {
  return StateModifyingFunctionMap();
}

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
