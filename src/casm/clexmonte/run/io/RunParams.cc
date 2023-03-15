#include "casm/clexmonte/run/io/RunParams.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
RunParams::RunParams(
    std::unique_ptr<state_generator_type> _state_generator,
    monte::RunManagerParams _run_manager_params,
    std::vector<sampling_figure_params_type> _sampling_fixture_params)
    : state_generator(std::move(_state_generator)),
      run_manager_params(std::move(_run_manager_params)),
      sampling_fixture_params(std::move(_sampling_fixture_params)) {}

}  // namespace clexmonte
}  // namespace CASM
