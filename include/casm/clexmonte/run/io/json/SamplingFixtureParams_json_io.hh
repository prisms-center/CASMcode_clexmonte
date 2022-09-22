#ifndef CASM_clexmonte_run_SamplingFixtureParams_json_io
#define CASM_clexmonte_run_SamplingFixtureParams_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

struct RunParams;

void parse(
    InputParser<monte::SamplingFixtureParams<config_type>> &parser,
    std::string label,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
    MethodParserMap<results_io_type> const &results_io_methods);

}  // namespace clexmonte
}  // namespace CASM

#endif
