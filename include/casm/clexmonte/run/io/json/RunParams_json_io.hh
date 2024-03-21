#ifndef CASM_clexmonte_run_RunParams_json_io
#define CASM_clexmonte_run_RunParams_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/global/filesystem.hh"
#include "casm/monte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

template <typename EngineType>
struct RunParams;

MethodParserMap<config_generator_type> standard_config_generator_methods(
    std::shared_ptr<system_type> const &system);

template <typename ConditionsType>
MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    StateModifyingFunctionMap const &modifying_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods,
    ConditionsType const *ptr = nullptr);

MethodParserMap<results_io_type> standard_results_io_methods();

template <typename EngineType, typename ConditionsType>
void parse(InputParser<RunParams<EngineType>> &parser,
           std::vector<fs::path> search_path,
           std::shared_ptr<EngineType> engine,
           monte::StateSamplingFunctionMap const &sampling_functions,
           monte::jsonStateSamplingFunctionMap const &json_sampling_functions,
           monte::ResultsAnalysisFunctionMap<config_type, statistics_type> const
               &analysis_functions,
           MethodParserMap<state_generator_type> const &state_generator_methods,
           MethodParserMap<results_io_type> const &results_io_methods,
           bool time_sampling_allowed, ConditionsType const *ptr = nullptr);

}  // namespace clexmonte
}  // namespace CASM

#endif
