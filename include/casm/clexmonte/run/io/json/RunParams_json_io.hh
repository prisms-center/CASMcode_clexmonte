#ifndef CASM_clexmonte_run_RunParams_json_io
#define CASM_clexmonte_run_RunParams_json_io

#include "casm/clexmonte/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

struct RunParams;

void parse(
    InputParser<RunParams> &parser, std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions);

}  // namespace clexmonte
}  // namespace CASM

#endif
