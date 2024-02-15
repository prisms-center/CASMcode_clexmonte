#ifndef CASM_clexmonte_run_functions
#define CASM_clexmonte_run_functions

#include "casm/casm_io/Log.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/io/json/ValueMap_json_io.hh"
#include "casm/monte/run_management/Results.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/SamplingFixture.hh"
#include "casm/monte/run_management/StateGenerator.hh"
#include "casm/monte/run_management/io/ResultsIO.hh"
#include "casm/monte/run_management/io/json/RunData_json_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief Perform a series of runs, according to a state_generator
template <typename CalculationType>
void run_series(
    CalculationType &calculation, state_generator_type &state_generator,
    monte::RunManagerParams const &run_manager_params,
    std::vector<sampling_fixture_params_type> const &sampling_fixture_params,
    std::vector<sampling_fixture_params_type> const &before_first_run =
        std::vector<sampling_fixture_params_type>({}),
    std::vector<sampling_fixture_params_type> const &before_each_run =
        std::vector<sampling_fixture_params_type>({}));

// --- Implementation ---

/// \brief Perform a series of runs, according to a state_generator
///
/// \param calculation A calculation instance, such as canonical::Canonical,
///     semigrand_canonical::SemiGrandCanonical, or kinetic::Kinetic.
/// \param state_generator A monte::StateGenerator, which produces a
///     a series of initial states
/// \param run_manager_params Parameters controlling the run manager
/// \param sampling_fixture_params Parameters controlling each
///     requested sampling fixture
///
/// Requires:
/// - std::shared_ptr<system_type> CalculationType::system: Shared ptr
///   with system info
/// - CalculationType::run(...): Method to run a single calculation, see
///   canonical::Canonical<EngineType>::run for an example
/// - bool CalculationType::update_species: For occupant tracking,
///   should be true for KMC, false otherwise
template <typename CalculationType>
void run_series(
    CalculationType &calculation, state_generator_type &state_generator,
    monte::RunManagerParams const &run_manager_params,
    std::vector<sampling_fixture_params_type> const &sampling_fixture_params,
    std::vector<sampling_fixture_params_type> const &before_first_run,
    std::vector<sampling_fixture_params_type> const &before_each_run) {
  typedef typename CalculationType::engine_type engine_type;
  std::shared_ptr<engine_type> engine =
      calculation.random_number_generator.engine;

  auto &log = CASM::log();
  log.begin("Monte Carlo calculation series");

  run_manager_type<engine_type> run_manager(run_manager_params, engine,
                                            sampling_fixture_params);

  // Final states are made available to the state generator which can use them
  // to determine the next state and enable restarts
  log.indent() << "Checking for completed runs..." << std::endl;
  run_manager.read_completed_runs();
  log.indent() << "Found " << run_manager.completed_runs.size() << std::endl
               << std::endl;

  // For all states generated, prepare input and run canonical Monte Carlo
  while (!state_generator.is_complete(run_manager.completed_runs)) {
    // Get initial state for the next calculation
    log.indent() << "Generating next state..." << std::endl;
    state_type state = state_generator.next_state(run_manager.completed_runs);
    log.indent() << qto_json(state.conditions) << std::endl;
    log.indent() << "Done" << std::endl;

    // Construct and initialize occupant tracking
    monte::Conversions const &convert =
        get_index_conversions(*calculation.system, state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(*calculation.system, state);
    monte::OccLocation occ_location(convert, occ_candidate_list,
                                    calculation.update_species);
    occ_location.initialize(get_occupation(state));

    // Optional, before first run:
    if (before_first_run.size() && run_manager.completed_runs.size() == 0) {
      run_manager_type<engine_type> tmp_run_manager(run_manager_params, engine,
                                                    before_first_run);
      // Run Monte Carlo at a single condition
      log.indent() << "Performing \"before-first-run\" run ..." << std::endl;
      calculation.run(state, occ_location, tmp_run_manager);
      log.indent() << "\"Before-first-run\" run: Done" << std::endl;
    }

    // Optional, before each run:
    if (before_each_run.size()) {
      run_manager_type<engine_type> tmp_run_manager(run_manager_params, engine,
                                                    before_each_run);
      // Run Monte Carlo at a single condition
      log.indent() << "Performing \"before-each-run\" run ..." << std::endl;
      calculation.run(state, occ_location, tmp_run_manager);
      log.indent() << "\"Before-each-run\" run: Done" << std::endl;
    }

    // Run Monte Carlo at a single condition
    log.indent() << "Performing Run " << run_manager.completed_runs.size() + 1
                 << "..." << std::endl;
    calculation.run(state, occ_location, run_manager);
    log.indent() << "Run " << run_manager.completed_runs.size() << " Done"
                 << std::endl;

    log.indent() << std::endl;
  }
  log.indent() << "Monte Carlo calculation series complete" << std::endl;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
