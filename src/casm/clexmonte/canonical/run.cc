#include "casm/clexmonte/canonical/run.hh"

#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/sampling_functions.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/methods/canonical.hh"
#include "casm/monte/results/Results.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/state/StateGenerator.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

void run(std::shared_ptr<system_type> const &system_data,
         state_generator_type &state_generator,
         monte::SamplingParams const &sampling_params,
         monte::CompletionCheckParams const &completion_check_params,
         monte::ResultsIO<config_type> &results_io,
         MTRand random_number_generator) {
  // Final states are used by the state generator to determine the next state
  std::vector<state_type> final_states;

  // Enable restarts: Check for a partically completed path
  final_states = results_io.read_final_states();

  // For all states generated, prepare input and run canonical Monte Carlo
  while (!state_generator.is_complete(final_states)) {
    // Get initial state for the next calculation
    state_type initial_state = state_generator.next_state(final_states);

    // Get supercell of initial state
    Eigen::Matrix3l const &transformation_matrix_to_super =
        get_transformation_matrix_to_super(initial_state.configuration);

    // Make supercell-specific formation energy clex calculator
    clexulator::ClusterExpansion &formation_energy_clex_calculator =
        get_formation_energy_clex(*system_data, initial_state);

    // Make state sampling functions, with current supercell-specific info
    monte::StateSamplingFunctionMap<config_type> sampling_functions =
        make_sampling_functions(system_data);

    // Prepare allowed occupation swaps
    monte::Conversions convert{*get_shared_prim(*system_data),
                               transformation_matrix_to_super};
    monte::OccCandidateList occ_candidate_list{convert};
    std::vector<monte::OccSwap> canonical_swaps =
        make_canonical_swaps(convert, occ_candidate_list);

    // Run Monte Carlo at a single condition
    results_type result = monte::canonical(
        initial_state, formation_energy_clex_calculator, convert,
        occ_candidate_list, canonical_swaps, random_number_generator,
        sampling_params, sampling_functions, completion_check_params);

    // Store final state for state generation input
    final_states.push_back(result.trajectory.back());

    // Write results for this condition
    Index run_index = final_states.size();
    results_io.write_initial_state(result.trajectory.front(), run_index);
    results_io.write_final_state(result.trajectory.back(), run_index);
    results_io.write_trajectory(result.trajectory, run_index);
    results_io.write_observations(result.sampled_data, run_index);
    results_io.write_completion_check_results(result.completion_check_results,
                                              run_index);
    results_io.write_summary(result.completion_check_results, run_index);
  }
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
