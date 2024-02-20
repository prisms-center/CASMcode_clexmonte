#ifndef CASM_clexmonte_run_ConfigGenerator_json_io
#define CASM_clexmonte_run_ConfigGenerator_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/monte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct ConfigGenerator from JSON
void parse(
    InputParser<config_generator_type> &parser,
    MethodParserMap<config_generator_type> const &config_generator_methods);

/// \brief Construct FixedConfigGenerator from JSON
void parse(InputParser<FixedConfigGenerator> &parser,
           std::shared_ptr<system_type> const &system);

}  // namespace clexmonte
}  // namespace CASM

#endif
