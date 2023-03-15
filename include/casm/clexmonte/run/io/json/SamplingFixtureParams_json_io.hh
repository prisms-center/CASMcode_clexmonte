#ifndef CASM_clexmonte_run_SamplingFixtureParams_json_io
#define CASM_clexmonte_run_SamplingFixtureParams_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

struct RunParams;

void parse(InputParser<sampling_figure_params_type> &parser, std::string label,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions,
           MethodParserMap<results_io_type> const &results_io_methods);

}  // namespace clexmonte
}  // namespace CASM

#endif
