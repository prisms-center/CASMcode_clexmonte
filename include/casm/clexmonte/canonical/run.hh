#ifndef CASM_clexmonte_canonical_run
#define CASM_clexmonte_canonical_run

#include "casm/clexmonte/canonical/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

void run(std::shared_ptr<system_type> const &system_data,
         state_generator_type &state_generator,
         monte::SamplingParams const &sampling_params,
         monte::CompletionCheckParams const &completion_check_params,
         monte::ResultsIO<config_type> &results_io,
         MTRand random_number_generator);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
