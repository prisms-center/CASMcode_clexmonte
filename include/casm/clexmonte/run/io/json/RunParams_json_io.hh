#ifndef CASM_clexmonte_run_RunParams_json_io
#define CASM_clexmonte_run_RunParams_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

struct RunParams;

MethodParserMap<config_generator_type> standard_config_generator_methods(
    std::shared_ptr<system_type> const &system);

MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods);

MethodParserMap<results_io_type> standard_results_io_methods(
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions);

void parse(
    InputParser<RunParams> &parser, std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
    MethodParserMap<state_generator_type> const &state_generator_methods,
    MethodParserMap<results_io_type> const &results_io_methods);

}  // namespace clexmonte
}  // namespace CASM

#endif
