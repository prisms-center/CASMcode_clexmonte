#include "casm/clexmonte/canonical/io/InputData.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Constructor
InputData::InputData(
    std::shared_ptr<system_type> _system_data,
    std::unique_ptr<state_generator_type> _state_generator,
    monte::StateSamplingFunctionMap<config_type> const &_sampling_functions,
    monte::SamplingParams const &_sampling_params,
    monte::CompletionCheckParams const &_completion_check_params,
    std::unique_ptr<results_io_type> _results_io,
    MTRand _random_number_generator)
    : system_data(_system_data),
      state_generator(std::move(_state_generator)),
      sampling_functions(_sampling_functions),
      sampling_params(_sampling_params),
      completion_check_params(_completion_check_params),
      results_io(std::move(_results_io)),
      random_number_generator(_random_number_generator) {}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
