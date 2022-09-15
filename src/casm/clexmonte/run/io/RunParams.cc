#include "casm/clexmonte/run/io/RunParams.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
RunParams::RunParams(
    monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &_analysis_functions,
    std::unique_ptr<state_generator_type> _state_generator,
    monte::SamplingParams const &_sampling_params,
    monte::CompletionCheckParams const &_completion_check_params,
    std::unique_ptr<results_io_type> _results_io, monte::MethodLog _method_log)
    : sampling_functions(_sampling_functions),
      analysis_functions(_analysis_functions),
      state_generator(std::move(_state_generator)),
      sampling_params(_sampling_params),
      completion_check_params(_completion_check_params),
      results_io(std::move(_results_io)),
      method_log(_method_log) {}

}  // namespace clexmonte
}  // namespace CASM
