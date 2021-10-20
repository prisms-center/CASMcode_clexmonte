#ifndef CASM_clexmonte_canonical_StateGenerator_json_io
#define CASM_clexmonte_canonical_StateGenerator_json_io

#include "casm/clexmonte/canonical/io/json/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct StateGenerator from JSON
void parse(
    InputParser<state_generator_type> &parser,
    std::shared_ptr<system_type> const &system_data,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    canonical_tag tag);

/// \brief Construct IncrementalConditionsStateGenerator from JSON
void parse(
    InputParser<incremental_state_generator_type> &parser,
    std::shared_ptr<system_type> const &system_data,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
