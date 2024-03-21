#ifndef CASM_clexmonte_run_StateGenerator_json_io
#define CASM_clexmonte_run_StateGenerator_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/monte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct StateGenerator from JSON
void parse(
    InputParser<state_generator_type> &parser,
    MethodParserMap<state_generator_type> const &state_generator_methods);

/// \brief Construct IncrementalConditionsStateGenerator from JSON
template <typename ConditionsType>
void parse(InputParser<IncrementalConditionsStateGenerator> &parser,
           std::shared_ptr<system_type> const &system,
           StateModifyingFunctionMap const &modifying_functions,
           MethodParserMap<config_generator_type> config_generator_methods,
           ConditionsType const *ptr = nullptr);

}  // namespace clexmonte
}  // namespace CASM

#endif
