#ifndef CASM_clexmonte_canonical_run
#define CASM_clexmonte_canonical_run

#include "casm/clexmonte/canonical/definitions.hh"
#include "casm/monte/MethodLog.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Run canonical Monte Carlo calculations
void run(std::shared_ptr<system_type> const &system_data,
         state_generator_type &state_generator,
         monte::StateSampler<config_type> &state_sampler,
         monte::CompletionCheck &completion_check,
         monte::ResultsIO<config_type> &results_io,
         MTRand &random_number_generator,
         monte::MethodLog method_log = monte::MethodLog());

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
