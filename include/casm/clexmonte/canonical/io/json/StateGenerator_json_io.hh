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
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    canonical_tag tag);

/// \brief Construct conditions (monte::VectorValueMap) from JSON
void parse(InputParser<monte::VectorValueMap> &parser,
           std::shared_ptr<system_type> const &system_data, canonical_tag tag);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
