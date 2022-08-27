#include "casm/clexmonte/canonical/io/InputData.hh"

#include "casm/clexmonte/canonical/run.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Constructor
InputData::InputData(
    std::shared_ptr<system_type> _system,
    std::unique_ptr<state_generator_type> _state_generator,
    monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
    monte::SamplingParams const &_sampling_params,
    monte::CompletionCheckParams const &_completion_check_params,
    std::unique_ptr<results_io_type> _results_io,
    MTRand _random_number_generator)
    : system(_system),
      state_generator(std::move(_state_generator)),
      sampling_functions(_sampling_functions),
      sampling_params(_sampling_params),
      completion_check_params(_completion_check_params),
      results_io(std::move(_results_io)),
      random_number_generator(_random_number_generator) {}

/// \brief Run canonical Monte Carlo calculations
void run(InputData &input_data) {
  // Create state sampler
  // - This object holds sampling functions, tracks the number of steps &
  //   passes, determines when samples are due, takes samples and holds data
  // - Use its member functions `increment_step` and `sample_data_if_due` to
  //   collect samples according to `sampling_params`
  monte::StateSampler<config_type> state_sampler(input_data.sampling_params,
                                                 input_data.sampling_functions);

  // Create CompletionCheck method
  // - This object checks for min/max cutoffs and automatic convergence
  monte::CompletionCheck completion_check(input_data.completion_check_params);

  run(input_data.system, *input_data.state_generator, state_sampler,
      completion_check, *input_data.results_io,
      input_data.random_number_generator);
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
