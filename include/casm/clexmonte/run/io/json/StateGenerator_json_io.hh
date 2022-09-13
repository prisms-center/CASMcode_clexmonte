#ifndef CASM_clexmonte_run_StateGenerator_json_io
#define CASM_clexmonte_run_StateGenerator_json_io

#include "casm/clexmonte/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct StateGenerator from JSON
void parse(
    InputParser<state_generator_type> &parser,
    std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions);

/// \brief Construct IncrementalConditionsStateGenerator from JSON
void parse(
    InputParser<monte::IncrementalConditionsStateGenerator<Configuration>>
        &parser,
    std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions);

}  // namespace clexmonte
}  // namespace CASM

#endif
