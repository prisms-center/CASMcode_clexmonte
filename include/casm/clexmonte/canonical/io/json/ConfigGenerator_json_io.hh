#ifndef CASM_clexmonte_ConfigGenerator_json_io
#define CASM_clexmonte_ConfigGenerator_json_io

#include "casm/clexmonte/canonical/io/json/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct ConfigGenerator from JSON
void parse(InputParser<config_generator_type> &parser,
           std::shared_ptr<system_type> const &system_data);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
