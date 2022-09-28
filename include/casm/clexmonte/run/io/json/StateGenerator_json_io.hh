#ifndef CASM_clexmonte_run_StateGenerator_json_io
#define CASM_clexmonte_run_StateGenerator_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct StateGenerator from JSON
void parse(
    InputParser<state_generator_type> &parser,
    MethodParserMap<state_generator_type> const &state_generator_methods);

/// \brief Construct IncrementalConditionsStateGenerator from JSON
void parse(
    InputParser<monte::IncrementalConditionsStateGenerator<Configuration>>
        &parser,
    std::shared_ptr<system_type> const &system,
    monte::StateModifyingFunctionMap<config_type> const &modifying_functions,
    MethodParserMap<config_generator_type> config_generator_methods);

}  // namespace clexmonte
}  // namespace CASM

#endif
