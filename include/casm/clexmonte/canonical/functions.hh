#ifndef CASM_clexmonte_canonical_functions
#define CASM_clexmonte_canonical_functions

#include "casm/clexmonte/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
monte::StateSamplingFunctionMap<Configuration> make_sampling_functions(
    std::shared_ptr<system_type> const &system);

/// \brief Construct functions that may be used to analyze Monte Carlo
///     calculation results
monte::ResultsAnalysisFunctionMap<Configuration> make_analysis_functions(
    std::shared_ptr<system_type> const &system);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
