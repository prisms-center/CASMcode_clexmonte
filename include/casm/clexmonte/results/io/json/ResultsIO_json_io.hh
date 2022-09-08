#ifndef CASM_clexmonte_results_ResultsIO_json_io
#define CASM_clexmonte_results_ResultsIO_json_io

#include "casm/monte/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct ResultsIO from JSON
template <typename ConfigType>
void parse(
    InputParser<monte::ResultsIO<ConfigType>> &parser,
    monte::StateSamplingFunctionMap<ConfigType> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<ConfigType> const &analysis_functions);

/// \brief Construct jsonResultsIO from JSON
template <typename ConfigType>
void parse(
    InputParser<monte::jsonResultsIO<ConfigType>> &parser,
    monte::StateSamplingFunctionMap<ConfigType> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<ConfigType> const &analysis_functions);

}  // namespace clexmonte
}  // namespace CASM

#endif
