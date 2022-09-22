#ifndef CASM_clexmonte_kinetic
#define CASM_clexmonte_kinetic

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/CompleteEventCalculator.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/EventStateCalculator.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Implements kinetic Monte Carlo calculations
template <typename EngineType>
struct Kinetic {
  explicit Kinetic(std::shared_ptr<system_type> _system,
                   std::shared_ptr<EngineType> _random_number_engine =
                       std::shared_ptr<EngineType>());

  /// System data
  std::shared_ptr<system_type> system;

  /// Random number generator
  monte::RandomNumberGenerator<EngineType> random_number_generator;

  /// Update species in monte::OccLocation tracker
  bool update_species = true;

  /// Current state
  monte::State<Configuration> const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_supercell;

  /// Occupant tracker
  monte::OccLocation const *occ_location;

  /// The current state's conditions in efficient-to-use form
  std::shared_ptr<clexmonte::Conditions> conditions;

  // TODO:
  // /// If true: rejection-free KMC, if false: rejection-KMC
  // bool rejection_free = true;

  // --- System specific ---

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  // --- Supercell & state specific ---

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Functions for calculating event states, one for each prim event.
  /// This is supercell-specific, even though it is one per prim event,
  /// because it depends on supercell-specific clexulators
  std::vector<clexmonte::EventStateCalculator> prim_event_calculators;

  /// Calculator for KMC event selection
  std::shared_ptr<clexmonte::CompleteEventCalculator> event_calculator;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           monte::RunManager<config_type> &run_manager);

  /// \brief Perform a series of runs, according to a state generator
  void run_series(state_generator_type &state_generator,
                  std::vector<monte::SamplingFixtureParams<config_type>> const
                      &sampling_fixture_params);

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  static monte::StateSamplingFunctionMap<Configuration>
  standard_sampling_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static monte::ResultsAnalysisFunctionMap<Configuration>
  standard_analysis_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);
};

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
