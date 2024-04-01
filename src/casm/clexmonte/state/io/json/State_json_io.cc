#include "casm/clexmonte/state/io/json/State_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/io/json/ValueMap_json_io.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {

/// \brief Write monte::State<clexmonte::Configuration> to JSON
jsonParser &to_json(monte::State<clexmonte::Configuration> const &state,
                    jsonParser &json, bool write_prim_basis) {
  to_json(state.configuration, json["configuration"], write_prim_basis);
  json["conditions"] = state.conditions;
  json["properties"] = state.properties;
  return json;
}

void parse(InputParser<monte::State<clexmonte::Configuration>> &parser,
           config::SupercellSet &supercells) {
  auto configuration_subparser =
      parser.subparse<clexmonte::Configuration>("configuration", supercells);
  auto conditions_subparser = parser.subparse<monte::ValueMap>("conditions");
  auto properties_subparser = parser.subparse<monte::ValueMap>("properties");

  if (parser.valid()) {
    parser.value = std::make_unique<monte::State<clexmonte::Configuration>>(
        *configuration_subparser->value, *conditions_subparser->value,
        *properties_subparser->value);
  }
}

/// \brief Read monte::State<clexmonte::Configuration> from JSON
monte::State<clexmonte::Configuration>
jsonConstructor<monte::State<clexmonte::Configuration>>::from_json(
    jsonParser const &json, config::SupercellSet &supercells) {
  InputParser<monte::State<clexmonte::Configuration>> parser{json, supercells};
  std::stringstream ss;
  ss << "Error reading monte::State<clexmonte::Configuration> from JSON input";
  report_and_throw_if_invalid(parser, CASM::log(),
                              std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Read monte::State<clexmonte::Configuration> from JSON
void from_json(monte::State<clexmonte::Configuration> &state,
               jsonParser const &json, config::SupercellSet &supercells) {
  state = jsonConstructor<monte::State<clexmonte::Configuration>>::from_json(
      json, supercells);
}

}  // namespace CASM
