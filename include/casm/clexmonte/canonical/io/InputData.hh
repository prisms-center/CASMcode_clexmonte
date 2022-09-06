#ifndef CASM_clexmonte_canonical_io_InputData
#define CASM_clexmonte_canonical_io_InputData

#include "casm/clexmonte/canonical/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/sampling/SamplingParams.hh"
#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Data structure for holding data parsed from the canonical Monte
///     Carlo input file
struct InputData {
  /// \brief Constructor
  InputData(
      std::shared_ptr<system_type> _system,
      std::unique_ptr<state_generator_type> _state_generator,
      monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
      monte::SamplingParams const &_sampling_params,
      monte::CompletionCheckParams const &_completion_check_params,
      std::unique_ptr<results_io_type> _results_io);

  /// System information:
  /// - prim
  /// - composition axes
  /// - formation energy clex data
  std::shared_ptr<system_type> system;

  /// State generator implementation
  std::unique_ptr<state_generator_type> state_generator;

  /// State sampling functions
  monte::StateSamplingFunctionMap<config_type> sampling_functions;

  /// Sampling parameters
  monte::SamplingParams sampling_params;

  /// Completion check params
  monte::CompletionCheckParams completion_check_params;

  /// Results I/O implementation
  std::unique_ptr<results_io_type> results_io;
};

/// \brief Run canonical Monte Carlo calculations
void run(InputData &input_data);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
