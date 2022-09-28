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

  // TODO:
  // /// If true: rejection-free KMC, if false: rejection-KMC
  // bool rejection_free = true;

  // --- System specific ---

  /// Specifies OccEvent index meaning / atom name indices
  std::shared_ptr<occ_events::OccSystem> event_system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  // --- Standard state specific ---

  /// Pointer to current state
  state_type const *state;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_supercell;

  /// Pointer to current occupant tracker
  monte::OccLocation const *occ_location;

  /// The current state's conditions in efficient-to-use form
  ///
  /// Note: This is shared with the calculators in `prim_event_calculators`
  std::shared_ptr<clexmonte::Conditions> conditions;

  /// When sampling, this will hold the atom name index for each column of the
  /// atom position matrices. Currently atom names only; does not distinguish
  /// atoms with different properties.
  std::vector<Index> atom_name_index_list;

  // TODO: need a better system for passing info to sampling functions

  /// When sampling, this will specify the current sampling fixture
  std::string sampling_fixture_label;

  /// When sampling, this will point to the current state sampler
  monte::StateSampler<clexmonte::config_type> const *state_sampler;

  /// When sampling, this will hold the current simulation time
  double time;

  /// When sampling, this will hold current atom positions
  Eigen::MatrixXd atom_positions_cart;

  /// When sampling, this will hold the last sampling time
  ///
  /// Notes:
  /// - Key = sampling fixture label
  /// - For the first sample, this will contain 0.0
  std::map<std::string, double> prev_time;

  /// When sampling, this will hold atom positions from the last sampling time
  ///
  /// Notes:
  /// - Key = sampling fixture label
  /// - For the first sample, this will contain the atom positions at the
  ///   start of the run.
  std::map<std::string, Eigen::MatrixXd> prev_atom_positions_cart;

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
  static monte::StateSamplingFunctionMap<config_type>
  standard_sampling_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  static monte::ResultsAnalysisFunctionMap<config_type>
  standard_analysis_functions(
      std::shared_ptr<Kinetic<EngineType>> const &calculation);

  /// \brief Construct functions that may be used to modify states
  static monte::StateModifyingFunctionMap<config_type>
  standard_modifying_functions(
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
