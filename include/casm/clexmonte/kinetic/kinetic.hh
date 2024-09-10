#ifndef CASM_clexmonte_kinetic
#define CASM_clexmonte_kinetic

#include <random>

#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Implements kinetic Monte Carlo calculations
template <typename EngineType>
struct Kinetic {
  typedef EngineType engine_type;

  explicit Kinetic(std::shared_ptr<system_type> _system,
                   std::vector<EventFilterGroup> _event_filters = {});

  /// System data
  std::shared_ptr<system_type> system;

  /// Event filters
  std::vector<EventFilterGroup> event_filters;

  /// Update species in monte::OccLocation tracker
  bool update_atoms = true;

  /// Method allows time-based sampling
  bool time_sampling_allowed = true;

  // TODO:
  // /// If true: rejection-free KMC, if false: rejection-KMC
  // bool rejection_free = true;

  /// \brief KMC event data and calculators
  std::shared_ptr<KineticEventData> event_data;

  // --- Standard state specific ---

  /// Pointer to current state
  state_type const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Pointer to current occupant tracker
  monte::OccLocation const *occ_location;

  /// The current state's conditions in efficient-to-use form
  ///
  /// Note: This is shared with the calculators in `prim_event_calculators`
  std::shared_ptr<clexmonte::Conditions> conditions;

  /// The current state's potential calculator, set
  ///    when the `run` method is called - for sampling function only
  std::shared_ptr<canonical::CanonicalPotential> potential;

  /// The current state's formation energy cluster expansion calculator, set
  ///    when the `run` method is called
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy;

  // --- Data used by kinetic sampling functions ---

  /// Data for sampling functions
  monte::KMCData<config_type, statistics_type, engine_type> kmc_data;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<EngineType> &run_manager);

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  static std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  static std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to modify states
  static StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);
};

/// \brief Construct a list of atom names corresponding to OccLocation atoms
std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system);

/// \brief Helper for making a conditions ValueMap for kinetic Monte
///     Carlo calculations
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions ValueMap for kinetic Monte
///     Carlo calculations
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Explicitly instantiated Kinetic calculator
typedef Kinetic<std::mt19937_64> Kinetic_mt19937_64;

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
