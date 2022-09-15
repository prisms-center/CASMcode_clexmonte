#ifndef CASM_clexmonte_canonical
#define CASM_clexmonte_canonical

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Implements potential for canonical Monte Carlo
class CanonicalPotential {
 public:
  CanonicalPotential(
      std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex);

  /// \brief Reset pointer to state currently being calculated
  void set(monte::State<Configuration> const *state);

  /// \brief Pointer to state currently being calculated
  monte::State<Configuration> const *get() const;

  /// \brief Calculate (extensive) cluster expansion value
  double extensive_value();

  /// \brief Calculate change in (extensive) cluster expansion value due to a
  ///     series of occupation changes
  double occ_delta_extensive_value(std::vector<Index> const &linear_site_index,
                                   std::vector<int> const &new_occ);

 private:
  /// State to use
  monte::State<Configuration> const *m_state;

  /// Formation energy cluster expansion calculator;
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
};

/// \brief Set potential calculator so it evaluates using `state`
void set(CanonicalPotential &potential,
         monte::State<Configuration> const &state);

typedef CanonicalPotential potential_type;

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Implements canonical Monte Carlo calculations
template <typename EngineType>
struct Canonical {
  explicit Canonical(std::shared_ptr<system_type> _system,
                     std::shared_ptr<EngineType> _random_number_engine =
                         std::shared_ptr<EngineType>());

  /// System data
  std::shared_ptr<system_type> system;

  /// Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator;

  /// \brief Make the Canonical potential calculator, for use with templated
  /// methods
  potential_type make_potential(state_type const &state) const;

  /// \brief Perform a single run, evolving current state
  ///
  /// Notes:
  /// - this->state and this->occ_location are evolved and end in modified
  /// states
  monte::Results<config_type> run(
      state_type &state, monte::OccLocation &occ_location,
      monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
      monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
      monte::SamplingParams const &sampling_params,
      monte::CompletionCheckParams const &completion_check_params,
      monte::MethodLog method_log = monte::MethodLog());
};

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
template <typename EngineType>
monte::StateSamplingFunctionMap<Configuration> standard_sampling_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation);

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
template <typename EngineType>
monte::ResultsAnalysisFunctionMap<Configuration> standard_analysis_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation);

// --- Implementation ---

template <typename EngineType>
Canonical<EngineType>::Canonical(
    std::shared_ptr<system_type> _system,
    std::shared_ptr<EngineType> _random_number_engine)
    : system(_system), random_number_generator(_random_number_engine) {
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing Canonical: no 'formation_energy' clex.");
  }
}

/// \brief Make the Canonical potential calculator, for use with templated
/// methods
template <typename EngineType>
potential_type Canonical<EngineType>::make_potential(
    state_type const &state) const {
  CanonicalPotential potential(get_clex(*system, state, "formation_energy"));
  set(potential, state);
  return potential;
}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
monte::Results<config_type> Canonical<EngineType>::run(
    state_type &state, monte::OccLocation &occ_location,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
    monte::SamplingParams const &sampling_params,
    monte::CompletionCheckParams const &completion_check_params,
    monte::MethodLog method_log) {
  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("mol_composition")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `mol_composition` conditions not set.");
  }

  CanonicalPotential potential = this->make_potential(state);

  /// \brief Construct swaps
  monte::Conversions const &convert = get_index_conversions(*system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*system, state);

  std::vector<monte::OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  std::vector<monte::OccSwap> grand_canonical_swaps =
      make_grand_canonical_swaps(convert, occ_candidate_list);

  monte::StateSampler<config_type> state_sampler(sampling_params,
                                                 sampling_functions);

  // Create CompletionCheck method
  // - This object checks for min/max cutoffs and automatic convergence
  monte::CompletionCheck completion_check(completion_check_params);

  // Enforce composition
  clexmonte::enforce_composition(
      get_occupation(state),
      state.conditions.vector_values.at("mol_composition"),
      get_composition_calculator(*system), grand_canonical_swaps, occ_location,
      random_number_generator);

  // Run Monte Carlo at a single condition
  typedef monte::RandomNumberGenerator<EngineType> generator_type;
  results_type result = monte::occupation_metropolis(
      state, occ_location, potential, canonical_swaps,
      monte::propose_canonical_event<generator_type>, random_number_generator,
      state_sampler, completion_check, analysis_functions, method_log);

  return result;
}

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param calculation Shared pointer to Canonical calculation, which
///     can be used by sampling functions to access system and calculation data
///     such as the prim, the cluster expansion, and the composition axes.
///
template <typename EngineType>
monte::StateSamplingFunctionMap<Configuration> standard_sampling_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  auto const &system = calculation->system;
  std::vector<monte::StateSamplingFunction<Configuration>> functions = {
      make_temperature_f(system),       make_mol_composition_f(system),
      make_param_composition_f(system), make_formation_energy_corr_f(system),
      make_formation_energy_f(system),  make_potential_energy_f(system)};

  monte::StateSamplingFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
};

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
template <typename EngineType>
monte::ResultsAnalysisFunctionMap<Configuration> standard_analysis_functions(
    std::shared_ptr<Canonical<EngineType>> const &calculation) {
  // TODO: heat_capacity, susc, etc.
  return monte::ResultsAnalysisFunctionMap<Configuration>();
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
