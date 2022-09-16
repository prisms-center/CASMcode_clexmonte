#ifndef CASM_clexmonte_semi_grand_canonical
#define CASM_clexmonte_semi_grand_canonical

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

/// \brief Implements potential for semi-grand canonical Monte Carlo
class SemiGrandCanonicalPotential {
 public:
  SemiGrandCanonicalPotential(std::shared_ptr<system_type> _system);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state);

  /// \brief Pointer to state currently being calculated
  state_type const *get() const;

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
  state_type const *m_state;

  /// Formation energy cluster expansion calculator, depends on current state
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;

  /// Number of unit cells, depends on current state
  double m_n_unitcells;

  /// Conditions, depends on current state
  std::shared_ptr<Conditions> m_conditions;

  /// Index conversions, depends on current state
  monte::Conversions const *m_convert;
};

/// \brief Set potential calculator so it evaluates using `state`
void set(SemiGrandCanonicalPotential &potential, state_type const &state);

typedef SemiGrandCanonicalPotential potential_type;

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

  /// \brief Make the SemiGrandCanonical potential calculator, for use
  ///     with templated methods
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
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static monte::ResultsAnalysisFunctionMap<Configuration>
  standard_analysis_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);
};

/// \brief Explicitly instantiated SemiGrandCanonical calculator
typedef SemiGrandCanonical<std::mt19937_64> SemiGrandCanonical_mt19937_64;

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
