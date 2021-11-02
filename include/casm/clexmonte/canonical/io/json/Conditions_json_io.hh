#ifndef CASM_clexmonte_canonical_Conditions_json_io
#define CASM_clexmonte_canonical_Conditions_json_io

#include "casm/clexmonte/canonical/io/json/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct conditions (monte::VectorValueMap) from JSON
void parse_conditions(InputParser<monte::VectorValueMap> &parser,
                      std::shared_ptr<system_type> const &system_data,
                      canonical_tag tag);

/// \brief Construct conditions increment (monte::VectorValueMap) from JSON
void parse_conditions_increment(InputParser<monte::VectorValueMap> &parser,
                                std::shared_ptr<system_type> const &system_data,
                                canonical_tag tag);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
