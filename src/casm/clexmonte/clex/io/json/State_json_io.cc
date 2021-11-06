#include "casm/clexmonte/clex/io/json/State_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/clex/io/json/Configuration_json_io.hh"
#include "casm/monte/state/State.hh"

namespace CASM {

/// \brief Write monte::State<clexmonte::Configuration> to JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values are written in the
///   basis in which they are provided
jsonParser &to_json(monte::State<clexmonte::Configuration> const &state,
                    jsonParser &json) {
  json["configuration"] = state.configuration;
  for (auto const &value : state.conditions) {
    to_json(value.second, json["conditions"][value.first],
            CASM::jsonParser::as_array());
  }
  for (auto const &value : state.properties) {
    to_json(value.second, json["properties"][value.first],
            CASM::jsonParser::as_array());
  }
  return json;
}

/// \brief Read monte::State<clexmonte::Configuration> from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
template <>
monte::State<clexmonte::Configuration>
from_json<monte::State<clexmonte::Configuration>>(jsonParser const &json) {
  ParentInputParser parser{json};

  std::unique_ptr<clexmonte::Configuration> configuration =
      parser.require<clexmonte::Configuration>("configuration");

  monte::VectorValueMap conditions;
  parser.optional(conditions, "conditions");

  monte::VectorValueMap properties;
  parser.optional(properties, "properties");

  auto &log = CASM::log();
  std::runtime_error error_if_invalid{
      "Error reading monte::State<clexmonte::Configuration> from JSON input"};
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  return monte::State<clexmonte::Configuration>(*configuration, conditions,
                                                properties);
}

/// \brief Read monte::State<clexmonte::Configuration> from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
void from_json(monte::State<clexmonte::Configuration> &state,
               jsonParser const &json) {
  state = from_json<monte::State<clexmonte::Configuration>>(json);
}

}  // namespace CASM
