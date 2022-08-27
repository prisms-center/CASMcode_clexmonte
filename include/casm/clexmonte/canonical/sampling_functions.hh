#ifndef CASM_clexmonte_canonical_sampling_functions
#define CASM_clexmonte_canonical_sampling_functions

#include "casm/clexmonte/canonical/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
monte::StateSamplingFunctionMap<Configuration> make_sampling_functions(
    std::shared_ptr<system_type> const &system,
    canonical_tag tag = canonical_tag());

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
