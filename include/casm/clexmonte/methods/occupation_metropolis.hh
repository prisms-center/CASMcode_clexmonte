/// An implementation of an occupation Metropolis Monte Carlo
/// main loop that makes use of the RunManager provided by
/// casm/monte/run_management to implement sampling
/// fixtures and results data structures and input/output
/// methods.

#ifndef CASM_clexmonte_methods_occupation_metropolis
#define CASM_clexmonte_methods_occupation_metropolis

#include <map>
#include <string>
#include <vector>

#include "casm/monte/Conversions.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/metropolis.hh"
#include "casm/monte/run_management/Results.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"
#include "casm/monte/sampling/SamplingParams.hh"

namespace CASM {
namespace clexmonte {

template <typename PotentialOccDeltaPerSupercellF,
          typename ProposeOccEventFuntionType,
          typename ApplyOccEventFuntionType, typename ConfigType,
          typename StatisticsType, typename EngineType>
void occupation_metropolis_v2(
    monte::State<ConfigType> &state, monte::OccLocation &occ_location,
    double temperature,
    PotentialOccDeltaPerSupercellF potential_occ_delta_per_supercell_f,
    ProposeOccEventFuntionType propose_event_f,
    ApplyOccEventFuntionType apply_event_f,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager);

// --- Implementation ---

/// \brief Run an occupation metropolis Monte Carlo calculation
///
/// \param state The state. Consists of both the initial
///     configuration and conditions. Conditions must include `temperature`
///     and any others required by `potential`.
/// \param occ_location An occupant location tracker, which enables efficient
///     event proposal. It must already be initialized with the input state.
/// \param temperature The temperature, in K.
/// \param potential_occ_delta_per_supercell_f A function, with signature
///     `double potential_occ_delta_per_supercell_f(OccEvent const &)`, which
///     calculates the change in potential energy due to a proposed event.
/// \param possible_swaps A vector of possible swap types,
///     indicated by the asymmetric unit index and occupant index of the
///     sites potentially being swapped. Typically constructed from
///     `make_canonical_swaps` which generates all possible canonical swaps, or
///     `make_semigrand_canonical_swaps` which generates all possible grand
///      canonical swaps. It can also be a subset to restrict which swaps are
///     allowed.
/// \param propose_event_f A function, with signature
///     `OccEvent const & propose_event_f(RandomNumberGenerator<EngineType>
///     &random_number_generator)`, which proposes an event.
/// \param apply_event_f A function, with signature
///     `void apply_event_f(OccEvent const &)`, which updates the state and
///     occ_location after an event is accepted.
/// \param run_manager Contains random number engine, sampling fixtures, and
///     after completion holds final results
///
template <typename PotentialOccDeltaPerSupercellF,
          typename ProposeOccEventFuntionType,
          typename ApplyOccEventFuntionType, typename ConfigType,
          typename StatisticsType, typename EngineType>
void occupation_metropolis_v2(
    monte::State<ConfigType> &state, monte::OccLocation &occ_location,
    double temperature,
    PotentialOccDeltaPerSupercellF potential_occ_delta_per_supercell_f,
    ProposeOccEventFuntionType propose_event_f,
    ApplyOccEventFuntionType apply_event_f,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager) {
  // # construct RandomNumberGenerator
  monte::RandomNumberGenerator<EngineType> random_number_generator(
      run_manager.engine);

  Index steps_per_pass = occ_location.mol_size();

  // Used within the main loop:
  double beta = 1.0 / (CASM::KB * temperature);
  double delta_potential_energy;

  // Main loop
  run_manager.initialize(steps_per_pass);
  run_manager.sample_data_by_count_if_due(state);
  while (!run_manager.is_complete()) {
    // Write run status, if due (check clocktime vs status log frequency, but
    // only after #samples or #count changes)
    run_manager.write_status_if_due();

    // Propose an event
    monte::OccEvent const &event = propose_event_f(random_number_generator);

    // Calculate change in potential energy (per_supercell) due to event
    delta_potential_energy = potential_occ_delta_per_supercell_f(event);

    // Accept or reject event
    bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                        random_number_generator);

    // Apply accepted event
    if (accept) {
      run_manager.increment_n_accept();
      apply_event_f(event);
    } else {
      run_manager.increment_n_reject();
    }

    // Increment count
    run_manager.increment_step();

    // Sample data, if a sample is due by count
    run_manager.sample_data_by_count_if_due(state);
  }

  run_manager.finalize(state);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
