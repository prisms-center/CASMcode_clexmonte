#include "casm/clexmonte/run/io/RunParams.hh"

#include "casm/clexmonte/canonical/run.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
RunParams::RunParams(
    monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &_analysis_functions,
    std::unique_ptr<state_generator_type> _state_generator,
    monte::SamplingParams const &_sampling_params,
    monte::CompletionCheckParams const &_completion_check_params,
    std::unique_ptr<results_io_type> _results_io)
    : sampling_functions(_sampling_functions),
      analysis_functions(_analysis_functions),
      state_generator(std::move(_state_generator)),
      sampling_params(_sampling_params),
      completion_check_params(_completion_check_params),
      results_io(std::move(_results_io)) {}

// /// \brief Run canonical Monte Carlo calculations
// void run(RunParams &run_params) {
//   // Create state sampler
//   // - This object holds sampling functions, tracks the number of steps &
//   //   passes, determines when samples are due, takes samples and holds data
//   // - Use its member functions `increment_step` and `sample_data_if_due` to
//   //   collect samples according to `sampling_params`
//   monte::StateSampler<config_type> state_sampler(run_params.sampling_params,
//                                                  run_params.sampling_functions);
//
//   // Create CompletionCheck method
//   // - This object checks for min/max cutoffs and automatic convergence
//   monte::CompletionCheck
//   completion_check(run_params.completion_check_params);
//
//   // Default construct an empty random number generator ptr
//   auto random_number_engine = std::make_shared<std::mt19937_64>();
//
//   run(run_params.system, *run_params.state_generator, state_sampler,
//       completion_check, run_params.analysis_functions,
//       *run_params.results_io, random_number_engine);
// }

}  // namespace clexmonte
}  // namespace CASM
