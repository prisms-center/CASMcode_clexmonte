#include "casm/clexmonte/canonical/io/json/Conditions_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/io/json/parse_composition.hh"
#include "casm/composition/CompositionConverter.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct conditions (monte::VectorValueMap) from JSON
///
/// Parses canonical Monte Carlo conditions from JSON. If successfully parsed,
/// `parser->value` will contain a monte::VectorValueMap with:
/// - "temperature": (size 1)
/// - "comp_n": (size = system components size)
///
/// Expected:
///
///   "temperature": number (required)
///     Temperature in K.
///
///   "comp_n": array of number or dict (optional)
///     Composition in number per primitive cell. May be:
///
///     - An array of number, specifying the number of each component per
///       primitive cell, interpreted using the order of the `"components"`
///       specified for composition axes. The size must match the number of
///       components.
///
///     - A dict, where the keys are the component names, and values are
///       the number of that component per primitive cell. All components in
///       the system must be included.
///
///   "comp_x": array of number or dict (optional)
///     Parametric composition, in terms of the chosen composition axes. Will
///     be converted to `"comp_n"`. May
///     be:
///
///     - An array of number, corresponding to
///       `[comp_x(a), comp_x(b), ...]`. The size must match the number of
///       composition axes.
///
///     - A dict, where the keys are the axes names ("a", "b", etc.), and
///       values are the corresponding parametric composition value.
///       All composition axes must be included.
///
///
void parse_conditions(
    InputParser<monte::VectorValueMap> &parser,
    composition::CompositionConverter const &composition_converter,
    canonical_tag tag) {
  double temperature;
  parser.require(temperature, "temperature");

  // note:
  // - this reads & validates "comp_n" or "comp_x" (exactly one allowed)
  //   as a composition value, using function `parse_composition_for_comp_n`
  // - if successful, `comp_n_subparser->value` contains an
  //   Eigen::VectorXd representing comp_n
  bool is_increment = false;
  auto comp_n_subparser = parser.parse_as_with<Eigen::VectorXd>(
      parse_composition_for_comp_n, composition_converter, is_increment);

  if (parser.valid()) {
    parser.value = std::make_unique<monte::VectorValueMap>();
    monte::VectorValueMap &conditions = *parser.value;
    conditions["temperature"] = to_VectorXd(temperature);
    conditions["comp_n"] = *comp_n_subparser->value;
  }
}

/// \brief Construct conditions increment (monte::VectorValueMap) from JSON
///
/// Parses canonical Monte Carlo conditions increments from JSON. If
/// successfully parsed, `parser->value` will contain a monte::VectorValueMap
/// with:
/// - "temperature": (size 1)
/// - "comp_n": (size = system components size)
///
/// The expected JSON format is the same as documented for `parse_conditions`,
/// but values are interpreted as increments.
void parse_conditions_increment(
    InputParser<monte::VectorValueMap> &parser,
    composition::CompositionConverter const &composition_converter,
    canonical_tag tag) {
  double temperature;
  parser.require(temperature, "temperature");

  // note:
  // - this reads & validates "comp_n" or "comp_x" (exactly one allowed)
  //   as a composition increment, using function `parse_composition_for_comp_n`
  // - if successful, `comp_n_subparser->value` contains an
  //   Eigen::VectorXd representing comp_n increment
  bool is_increment = true;
  auto comp_n_subparser = parser.parse_as_with<Eigen::VectorXd>(
      parse_composition_for_comp_n, composition_converter, is_increment);

  if (parser.valid()) {
    parser.value = std::make_unique<monte::VectorValueMap>();
    monte::VectorValueMap &conditions = *parser.value;
    conditions["temperature"] = to_VectorXd(temperature);
    conditions["comp_n"] = *comp_n_subparser->value;
  }
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
