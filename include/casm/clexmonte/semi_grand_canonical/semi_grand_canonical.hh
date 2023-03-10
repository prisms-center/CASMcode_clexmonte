#ifndef CASM_clexmonte_semi_grand_canonical
#define CASM_clexmonte_semi_grand_canonical

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical_events.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

/// \brief Helper for making a conditions ValueMap for semi-grand
///     canonical Monte Carlo calculations
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> param_chem_pot);

/// \brief Helper for making a conditions increment ValueMap for
///     semi-grand canonical Monte Carlo calculations
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> param_chem_pot);

/// \brief Implements semi-grand canonical Monte Carlo calculations
template <typename EngineType>
struct SemiGrandCanonical {
  explicit SemiGrandCanonical(
      std::shared_ptr<system_type> _system,
      std::shared_ptr<EngineType> _random_number_engine =
          std::shared_ptr<EngineType>());

  /// System data
  std::shared_ptr<system_type> system;

  /// Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator;

  /// Update species in monte::OccLocation tracker?
  bool update_species = false;

  /// Current state
  monte::State<Configuration> const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Occupant tracker
  monte::OccLocation const *occ_location;

  /// Pointer to current state sampler
  monte::StateSampler<Configuration> const *state_sampler;

  /// The current state's conditions in efficient-to-use form
  std::shared_ptr<clexmonte::Conditions> conditions;

  /// Data for N-fold way implementation
  std::shared_ptr<SemiGrandCanonicalEventData> event_data;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           monte::RunManager<config_type> &run_manager);

  /// \brief Perform a series of runs, according to a state generator
  void run_series(state_generator_type &state_generator,
                  std::vector<monte::SamplingFixtureParams<config_type>> const
                      &sampling_fixture_params);

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  static monte::StateSamplingFunctionMap<Configuration>
  standard_sampling_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static monte::ResultsAnalysisFunctionMap<Configuration>
  standard_analysis_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to modify states
  static monte::StateModifyingFunctionMap<config_type>
  standard_modifying_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);
};

/// \brief Explicitly instantiated SemiGrandCanonical calculator
typedef SemiGrandCanonical<std::mt19937_64> SemiGrandCanonical_mt19937_64;

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
