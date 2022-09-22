#ifndef CASM_clexmonte_run_io_RunParams
#define CASM_clexmonte_run_io_RunParams

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/SamplingFixture.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/results/ResultsAnalysisFunction.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/sampling/SamplingParams.hh"
#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

/// \brief Data structure for holding data parsed from the canonical Monte
///     Carlo input file
struct RunParams {
  /// \brief Constructor
  RunParams(std::unique_ptr<state_generator_type> _state_generator,
            std::vector<monte::SamplingFixtureParams<config_type>>
                _sampling_fixture_params);

  /// State generator implementation
  std::unique_ptr<state_generator_type> state_generator;

  /// Parameters for 0 or more sampling fixtures
  std::vector<monte::SamplingFixtureParams<config_type>>
      sampling_fixture_params;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
