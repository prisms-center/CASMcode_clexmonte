#ifndef CASM_clexmonte_run_RunParams_json_io
#define CASM_clexmonte_run_RunParams_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

template <typename EngineType>
struct RunParams;

MethodParserMap<config_generator_type> standard_config_generator_methods(
    std::shared_ptr<system_type> const &system);

MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    std::map<std::string, state_modifying_function_type> const
        &modifying_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods);

MethodParserMap<results_io_type> standard_results_io_methods(
    std::map<std::string, state_sampling_function_type> const
        &sampling_functions,
    std::map<std::string, results_analysis_function_type> const
        &analysis_functions);

template <typename EngineType>
void parse(InputParser<RunParams<EngineType>> &parser,
           std::shared_ptr<EngineType> engine,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions,
           MethodParserMap<state_generator_type> const &state_generator_methods,
           MethodParserMap<results_io_type> const &results_io_methods,
           bool time_sampling_allowed);

}  // namespace clexmonte
}  // namespace CASM

#endif
