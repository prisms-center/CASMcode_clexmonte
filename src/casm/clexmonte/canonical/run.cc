#include "casm/clexmonte/canonical/run.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/CanonicalPotential.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/clex/io/json/State_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/enforce_composition.hh"
#include "casm/clexmonte/system/sampling_functions.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/methods/occupation_metropolis.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Run canonical Monte Carlo calculations
///
/// Required interfaces:
/// - get_formation_energy_clex(system_type, state_type)
/// - get_shared_prim(system_type)
/// - get_composition_calculator(system_type)
/// - get_transformation_matrix_to_super(config_type)
/// - get_occupation(config_type)
///
/// Required state conditions:
/// - `temperature`: size=1
///   The temperature in K.
/// - `comp_n`: size=n_components
///   Size must match `system_data->composition_converter.components().size()`.
///
/// State properties that are set:
/// - `potential_energy`: size=1
///   The intensive potential energy (eV / unit cell).
///
void run(std::shared_ptr<system_type> const &system_data,
         state_generator_type &state_generator,
         monte::StateSampler<config_type> &state_sampler,
         monte::CompletionCheck &completion_check,
         monte::ResultsIO<config_type> &results_io,
         MTRand &random_number_generator) {
  auto &log = CASM::log();
  log.begin("Cluster expansion canonical Monte Carlo");

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
    log.indent() << as_flattest_json(initial_state.conditions) << std::endl;
    log.indent() << "Done" << std::endl;

    // Make supercell-specific potential calculator
    CanonicalPotential potential(
        get_formation_energy_clex(*system_data, initial_state));

    // Prepare canonical swaps -- currently all allowed, but could be selected
    monte::Conversions convert =
        get_index_conversions(*system_data, initial_state);
    monte::OccCandidateList const &occ_candidate_list =
        get_occ_candidate_list(*system_data, initial_state);
    std::vector<monte::OccSwap> canonical_swaps =
        make_canonical_swaps(convert, occ_candidate_list);
    std::vector<monte::OccSwap> grand_canonical_swaps =
        make_grand_canonical_swaps(convert, occ_candidate_list);

    log.indent() << "Enforcing composition..." << std::endl;
    enforce_composition(get_occupation(initial_state.configuration),
                        initial_state.conditions.at("mol_composition"),
                        get_composition_calculator(*system_data), convert,
                        grand_canonical_swaps, random_number_generator);
    log.indent() << "Done" << std::endl;

    // Run Monte Carlo at a single condition
    log.indent() << "Beginning run " << final_states.size() + 1 << std::endl;
    results_type result = monte::occupation_metropolis(
        initial_state, potential, convert, canonical_swaps,
        monte::propose_canonical_event, random_number_generator, state_sampler,
        completion_check);
    log.indent() << "Run complete" << std::endl;

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
