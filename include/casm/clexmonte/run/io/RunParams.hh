#ifndef CASM_clexmonte_run_io_RunParams
#define CASM_clexmonte_run_io_RunParams

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/RunManager.hh"
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
template <typename EngineType>
struct RunParams {
  /// \brief Constructor
  RunParams(std::shared_ptr<EngineType> _engine,
            std::unique_ptr<state_generator_type> _state_generator,
            monte::RunManagerParams _run_manager_params,
            std::vector<sampling_fixture_params_type> _sampling_fixture_params,
            std::vector<sampling_fixture_params_type> _before_first_run =
                std::vector<sampling_fixture_params_type>({}),
            std::vector<sampling_fixture_params_type> _before_each_run =
                std::vector<sampling_fixture_params_type>({}));

  /// Random number generator engine
  std::shared_ptr<EngineType> engine;

  /// State generator implementation
  std::unique_ptr<state_generator_type> state_generator;

  /// Run manager parameters
  monte::RunManagerParams run_manager_params;

  /// Parameters for 0 or more sampling fixtures
  std::vector<sampling_fixture_params_type> sampling_fixture_params;

  std::vector<sampling_fixture_params_type> before_first_run;

  std::vector<sampling_fixture_params_type> before_each_run;
};

/// --- template implementation ---

/// \brief Constructor
template <typename EngineType>
RunParams<EngineType>::RunParams(
    std::shared_ptr<EngineType> _engine,
    std::unique_ptr<state_generator_type> _state_generator,
    monte::RunManagerParams _run_manager_params,
    std::vector<sampling_fixture_params_type> _sampling_fixture_params,
    std::vector<sampling_fixture_params_type> _before_first_run,
    std::vector<sampling_fixture_params_type> _before_each_run)
    : engine(_engine),
      state_generator(std::move(_state_generator)),
      run_manager_params(std::move(_run_manager_params)),
      sampling_fixture_params(std::move(_sampling_fixture_params)),
      before_first_run(std::move(_before_first_run)),
      before_each_run(std::move(_before_each_run)) {}

}  // namespace clexmonte
}  // namespace CASM

#endif
