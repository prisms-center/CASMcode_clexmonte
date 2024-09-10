#ifndef CASM_clexmonte_semigrand_canonical
#define CASM_clexmonte_semigrand_canonical

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexmonte/semigrand_canonical/potential.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

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
  typedef EngineType engine_type;
  typedef SemiGrandCanonicalConditions conditions_type;

  explicit SemiGrandCanonical(std::shared_ptr<system_type> _system);

  /// System data
  std::shared_ptr<system_type> system;

  /// Update species in monte::OccLocation tracker?
  bool update_atoms = false;

  /// Method does not allow time-based sampling
  bool time_sampling_allowed = false;

  /// Current state
  state_type const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Occupant tracker
  monte::OccLocation const *occ_location;

  /// The current state's conditions in efficient-to-use form, set
  ///    when the `run` method is called
  std::shared_ptr<SemiGrandCanonicalConditions> conditions;

  /// The current state's potential calculator, set
  ///    when the `run` method is called
  std::shared_ptr<SemiGrandCanonicalPotential> potential;

  /// The current state's formation energy cluster expansion calculator, set
  ///    when the `run` method is called
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy;

  /// State sampling functions
  std::map<std::string, state_sampling_function_type> sampling_functions;

  /// JSON State sampling functions
  std::map<std::string, json_state_sampling_function_type>
      json_sampling_functions;

  /// Results analysis functions
  std::map<std::string, results_analysis_function_type> analysis_functions;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<EngineType> &run_manager);

  /// \brief Construct functions that may be used to sample various quantities
  /// of
  ///     the Monte Carlo calculation as it runs
  static std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  static std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to modify states
  static StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<SemiGrandCanonical<EngineType>> const &calculation);
};

/// \brief Explicitly instantiated SemiGrandCanonical calculator
typedef SemiGrandCanonical<std::mt19937_64> SemiGrandCanonical_mt19937_64;

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
