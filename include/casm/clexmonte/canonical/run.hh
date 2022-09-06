#ifndef CASM_clexmonte_canonical_run
#define CASM_clexmonte_canonical_run

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/CanonicalPotential.hh"
#include "casm/clexmonte/canonical/definitions.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"
#include "casm/monte/state/io/json/ValueMap_json_io.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Run canonical Monte Carlo calculations
template <typename EngineType>
void run(std::shared_ptr<system_type> const &system,
         state_generator_type &state_generator,
         monte::StateSampler<config_type> &state_sampler,
         monte::CompletionCheck &completion_check,
         monte::ResultsIO<config_type> &results_io,
         std::shared_ptr<EngineType> &random_number_engine =
             std::shared_ptr<EngineType>(),
         monte::MethodLog method_log = monte::MethodLog());

// --- Implementation ---

/// \brief Run canonical Monte Carlo calculations
///
/// \param random_number_engine A random number generating engine. Will be
///     default constructed and seeded by std::random_device if empty.
///
/// Required interfaces:
/// - get_clex(system_type, state_type, key)
/// - get_shared_prim(system_type)
/// - get_composition_calculator(system_type)
/// - get_transformation_matrix_to_super(config_type)
/// - get_occupation(config_type)
///
/// Required state conditions:
/// - `temperature`: size=1
///   The temperature in K.
/// - `mol_composition`: size=n_components
///   Size must match `system->composition_converter.components().size()`.
///
/// State properties that are set:
/// - `potential_energy`: size=1
///   The intensive potential energy (eV / unit cell).
///
template <typename EngineType>
void run(std::shared_ptr<system_type> const &system,
         state_generator_type &state_generator,
         monte::StateSampler<config_type> &state_sampler,
         monte::CompletionCheck &completion_check,
         monte::ResultsIO<config_type> &results_io,
         std::shared_ptr<EngineType> &random_number_engine,
         monte::MethodLog method_log) {
  auto &log = CASM::log();
  log.begin("Cluster expansion canonical Monte Carlo");

  typedef monte::RandomNumberGenerator<EngineType> generator_type;
  generator_type random_number_generator(random_number_engine);

  // Final states are made available to the state generator which can use them
  // to determine the next state
  std::vector<state_type> final_states;

  // Enable restarts: Check for a partically completed path
  log.indent() << "Checking for finished runs..." << std::endl;
  final_states = results_io.read_final_states();
  log.indent() << "Found " << final_states.size() << std::endl << std::endl;

  // For all states generated, prepare input and run canonical Monte Carlo
  while (!state_generator.is_complete(final_states)) {
    log.indent() << "Generating next initial state..." << std::endl;
    // Get initial state for the next calculation
    state_type initial_state = state_generator.next_state(final_states);
    log.indent() << to_json(initial_state.conditions) << std::endl;
    log.indent() << "Done" << std::endl;

    // Make supercell-specific potential calculator
    CanonicalPotential potential(
        get_clex(*system, initial_state, "formation_energy"));

    // Prepare canonical swaps -- currently all allowed, but could be selected
    monte::Conversions convert = get_index_conversions(*system, initial_state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(*system, initial_state);
    monte::OccLocation occ_location(convert, occ_candidate_list);
    std::vector<monte::OccSwap> canonical_swaps =
        make_canonical_swaps(convert, occ_candidate_list);
    std::vector<monte::OccSwap> grand_canonical_swaps =
        make_grand_canonical_swaps(convert, occ_candidate_list);

    log.indent() << "Enforcing composition..." << std::endl;
    enforce_composition(
        get_occupation(initial_state),
        initial_state.conditions.vector_values.at("mol_composition"),
        get_composition_calculator(*system), grand_canonical_swaps,
        occ_location, random_number_generator);
    log.indent() << "Done" << std::endl;

    // Run Monte Carlo at a single condition
    log.indent() << "Performing Run " << final_states.size() + 1 << "..."
                 << std::endl;
    results_type result = monte::occupation_metropolis(
        initial_state, potential, convert, canonical_swaps,
        monte::propose_canonical_event<generator_type>, random_number_generator,
        state_sampler, completion_check, method_log);
    log.indent() << "Run " << final_states.size() + 1 << " Done" << std::endl;

    // Store final state for state generation input
    final_states.push_back(*result.final_state);

    // Write results for this condition
    Index run_index = final_states.size();
    results_io.write(result, run_index);

    log.indent() << std::endl;
  }
  log.indent() << "Canonical Monte Carlo Done" << std::endl;
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
