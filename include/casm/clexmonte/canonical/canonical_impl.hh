#ifndef CASM_clexmonte_canonical_impl
#define CASM_clexmonte_canonical_impl

#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/modifying_functions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/run_management/Results.hh"
#include "casm/monte/run_management/State.hh"

// temporary
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/monte/events/io/OccCandidate_json_io.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

template <typename EngineType>
Canonical<EngineType>::Canonical(std::shared_ptr<system_type> _system)
    : system(_system),
      state(nullptr),
      transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)),
      occ_location(nullptr) {
  if (!is_clex_data(*this->system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing Canonical: no 'formation_energy' clex.");
  }
}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void Canonical<EngineType>::run(state_type &state,
                                monte::OccLocation &occ_location,
                                run_manager_type<EngineType> &run_manager) {
  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("mol_composition")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `mol_composition` conditions not set.");
  }

  this->state = &state;
  this->transformation_matrix_to_super =
      get_transformation_matrix_to_super(state);
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);

  // Make potential calculator
  this->potential = std::make_shared<CanonicalPotential>(this->system);
  this->potential->set(this->state, this->conditions);
  this->formation_energy = this->potential->formation_energy();

  /// \brief Get swaps
  std::vector<monte::OccSwap> const &canonical_swaps =
      get_canonical_swaps(*this->system);

  std::vector<monte::OccSwap> const &semigrand_canonical_swaps =
      get_semigrand_canonical_swaps(*this->system);

  // Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator(
      run_manager.engine);

  // Enforce composition
  clexmonte::enforce_composition(
      get_occupation(state),
      state.conditions.vector_values.at("mol_composition"),
      get_composition_calculator(*this->system), semigrand_canonical_swaps,
      occ_location, random_number_generator);

  // Run Monte Carlo at a single condition
  typedef monte::RandomNumberGenerator<EngineType> generator_type;
  monte::occupation_metropolis(state, occ_location, *this->potential,
                               canonical_swaps,
                               monte::propose_canonical_event<generator_type>,
                               random_number_generator, run_manager);
}

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param calculation Shared pointer to Canonical calculation, which
///     can be used by sampling functions to access system and calculation data
///     such as the prim, the cluster expansion, and the composition axes.
///
template <typename EngineType>
std::map<std::string, state_sampling_function_type>
Canonical<EngineType>::standard_sampling_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<state_sampling_function_type> functions = {
      make_temperature_f(calculation),
      make_mol_composition_f(calculation),
      make_param_composition_f(calculation),
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
Canonical<EngineType>::standard_json_sampling_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
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
Canonical<EngineType>::standard_analysis_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<results_analysis_function_type> functions = {
      make_heat_capacity_f(calculation)};

  std::map<std::string, results_analysis_function_type> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

/// \brief Construct functions that may be used to modify states
template <typename EngineType>
StateModifyingFunctionMap Canonical<EngineType>::standard_modifying_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  std::vector<StateModifyingFunction> functions = {
      make_set_mol_composition_f(calculation)};

  StateModifyingFunctionMap function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
