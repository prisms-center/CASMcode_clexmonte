#include "casm/clexmonte/run/io/RunParams.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
RunParams::RunParams(std::unique_ptr<state_generator_type> _state_generator,
                     std::vector<monte::SamplingFixtureParams<config_type>>
                         _sampling_fixture_params)
    : state_generator(std::move(_state_generator)),
      sampling_fixture_params(std::move(_sampling_fixture_params)) {}

}  // namespace clexmonte
}  // namespace CASM
