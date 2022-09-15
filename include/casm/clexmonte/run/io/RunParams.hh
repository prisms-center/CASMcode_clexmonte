#ifndef CASM_clexmonte_run_io_RunParams
#define CASM_clexmonte_run_io_RunParams

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/MethodLog.hh"
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
  RunParams(
      monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
      monte::ResultsAnalysisFunctionMap<config_type> const &_analysis_functions,
      std::unique_ptr<state_generator_type> _state_generator,
      monte::SamplingParams const &_sampling_params,
      monte::CompletionCheckParams const &_completion_check_params,
      std::unique_ptr<results_io_type> _results_io,
      monte::MethodLog _method_log);

  /// State sampling functions
  monte::StateSamplingFunctionMap<config_type> sampling_functions;

  /// Results analysis functions
  monte::ResultsAnalysisFunctionMap<config_type> analysis_functions;

  /// State generator implementation
  std::unique_ptr<state_generator_type> state_generator;

  /// Sampling parameters
  monte::SamplingParams sampling_params;

  /// Completion check params
  monte::CompletionCheckParams completion_check_params;

  /// Results I/O implementation
  std::unique_ptr<results_io_type> results_io;

  /// Logging
  monte::MethodLog method_log;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
