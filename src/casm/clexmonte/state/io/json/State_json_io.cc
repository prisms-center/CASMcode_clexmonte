#include "casm/clexmonte/state/io/json/State_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/monte/io/json/ValueMap_json_io.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {

/// \brief Write monte::State<clexmonte::Configuration> to JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values are written in the
///   basis in which they are provided
jsonParser &to_json(monte::State<clexmonte::Configuration> const &state,
                    jsonParser &json) {
  json["configuration"] = state.configuration;
  json["conditions"] = state.conditions;
  json["properties"] = state.properties;
  return json;
}

void parse(InputParser<monte::State<clexmonte::Configuration>> &parser) {
  std::unique_ptr<clexmonte::Configuration> configuration =
      parser.require<clexmonte::Configuration>("configuration");

  auto conditions_subparser = parser.subparse<monte::ValueMap>("conditions");
  auto properties_subparser = parser.subparse<monte::ValueMap>("properties");

  if (parser.valid()) {
    parser.value = std::make_unique<monte::State<clexmonte::Configuration>>(
        *configuration, *conditions_subparser->value,
        *properties_subparser->value);
  }
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
  InputParser<monte::State<clexmonte::Configuration>> parser{json};
  std::stringstream ss;
  ss << "Error reading monte::State<clexmonte::Configuration> from JSON input";
  report_and_throw_if_invalid(parser, CASM::log(),
                              std::runtime_error{ss.str()});
  return *parser.value;
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
