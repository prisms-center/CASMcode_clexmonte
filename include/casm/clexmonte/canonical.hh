#ifndef CASM_clexmonte_canonical
#define CASM_clexmonte_canonical

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Implements potential for canonical Monte Carlo
class CanonicalPotential {
 public:
  CanonicalPotential(std::shared_ptr<system_type> _system);

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
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// State to use
  monte::State<Configuration> const *m_state;

  /// Conditions, depends on current state
  std::shared_ptr<Conditions> m_conditions;

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

  /// Update species in monte::OccLocation tracker?
  bool update_species = false;

  /// \brief Make the Canonical potential calculator, for use with templated
  /// methods
  potential_type make_potential(state_type const &state) const;

  /// \brief Perform a single run, evolving current state
  monte::Results<config_type> run(
      state_type &state, monte::OccLocation &occ_location,
      monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
      monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
      monte::SamplingParams const &sampling_params,
      monte::CompletionCheckParams const &completion_check_params,
      monte::MethodLog method_log = monte::MethodLog());

  /// \brief Perform a series of runs, according to a state generator
  void run_series(
      monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
      monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
      monte::SamplingParams const &sampling_params,
      monte::CompletionCheckParams &completion_check_params,
      state_generator_type &state_generator, results_io_type &results_io,
      monte::MethodLog method_log = monte::MethodLog());

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  static monte::StateSamplingFunctionMap<Configuration>
  standard_sampling_functions(
      std::shared_ptr<Canonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static monte::ResultsAnalysisFunctionMap<Configuration>
  standard_analysis_functions(
      std::shared_ptr<Canonical<EngineType>> const &calculation);
};

/// \brief Explicitly instantiated Canonical calculator
typedef Canonical<std::mt19937_64> Canonical_mt19937_64;

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
