#ifndef CASM_clexmonte_semigrand_canonical_impl
#define CASM_clexmonte_semigrand_canonical_impl

#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/run/analysis_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/semigrand_canonical/semigrand_canonical.hh"
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
#include "casm/monte/run_management/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

template <typename EngineType>
SemiGrandCanonical<EngineType>::SemiGrandCanonical(
    std::shared_ptr<system_type> _system,
    std::shared_ptr<EngineType> _random_number_engine)
    : system(_system),
      random_number_generator(_random_number_engine),
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
  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("param_chem_pot")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `param_chem_pot` conditions not set.");
  }

  this->state = &state;
  this->transformation_matrix_to_super =
      get_transformation_matrix_to_super(state);
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);

  auto potential = std::make_shared<SemiGrandCanonicalPotential>(this->system);
  potential->set(this->state, this->conditions);

  /// \brief Construct swaps
  monte::Conversions const &convert =
      get_index_conversions(*this->system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*this->system, state);

  std::vector<monte::OccSwap> semigrand_canonical_swaps =
      make_semigrand_canonical_swaps(convert, occ_candidate_list);

  // Run Monte Carlo at a single condition
  typedef monte::RandomNumberGenerator<EngineType> generator_type;
  monte::occupation_metropolis(
      state, occ_location, *potential, semigrand_canonical_swaps,
      monte::propose_semigrand_canonical_event<generator_type>,
      random_number_generator, run_manager);
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
};

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
std::map<std::string, state_modifying_function_type>
SemiGrandCanonical<EngineType>::standard_modifying_functions(
    std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation) {
  return std::map<std::string, state_modifying_function_type>();
}

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
