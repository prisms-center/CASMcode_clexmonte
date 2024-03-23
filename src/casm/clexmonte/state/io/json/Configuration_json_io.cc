#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"

namespace CASM {

/// \defgroup clexmonte_Configuration_JSON clexmonte::Configuration JSON format
///
/// ### JSON Attributes List
///
/// | Name | Description | Format |
/// |-|-|-|
/// | `transformation_matrix_to_supercell` | Supercell defining transformation
/// matrix | 2d array of int | | `dof` | Configuration DoF values |
/// clexulator::ConfigDoFValues |
///
/// ### JSON Attributes Description
///
/// #### clexmonte::Configuration JSON object
///
/// - `transformation_matrix_to_supercell`: 2d array of int (required,
///       shape=(3,3))
///   - A 3x3 integer transformation matrix, \f$T\f$, that gives the supercell
///     lattice vectors, \f$L^{scel}\f$, in terms of the prim lattice vectors,
///     \f$L^{prim}\f$, according to:
///     \f[
///         L^{scel} = L^{prim} T,
///     \f]
///     where \f$L^{scel}\f$ and \f$L^{prim}\f$ are 3x3 matrices whose columns
///     are the lattice vectors.
/// - `dof`: clexulator::ConfigDoFValues
///   - Configuration DoF values. Depending on context may be expressed in the
///     standard or prim basis.
/// - `basis`: str
////  - One of "prim" or "standard"
///
/// ### Examples
///
/// #### Example 1) Configuration with occupation, displacement, and strain DoF
/// ```
/// {
///   "transformation_matrix_to_supercell": [],
///   "dof": {}
/// }
/// ```

/// \brief Write clexmonte::Configuration to JSON
///
/// Notes:
/// - This assumes the DoF values basis is "prim", and writes the DoF values
///   without any conversion
///
/// \ingroup clexmonte_Configuration_JSON
jsonParser &to_json(clexmonte::Configuration const &configuration,
                    jsonParser &json) {
  json["transformation_matrix_to_supercell"] =
      configuration.transformation_matrix_to_super;
  json["dof"] = configuration.dof_values;
  json["basis"] = "prim";
  return json;
}

void parse(InputParser<clexmonte::Configuration> &parser) {
  if (parser.self.contains("basis")) {
    if (!parser.self["basis"].is_string() ||
        parser.self["basis"].get<std::string>() != "prim") {
      parser.error.insert(
          "Error reading clexmonte::Configuration: "
          "The \"basis\" must be \"prim\".");
    }
  } else {
    parser.error.insert(
        "Error reading clexmonte::Configuration: "
        "Missing \"basis\" value, which must be "
        "\"prim\".");
  }
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  auto dof_values_subparser =
      parser.subparse<clexulator::ConfigDoFValues>("dof");

  if (parser.valid()) {
    parser.value = std::make_unique<clexmonte::Configuration>(
        T, *dof_values_subparser->value);
  }
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - The "basis" must be "prim"
/// - This does not check the validity of the DoF values dimensions
////
/// \ingroup clexmonte_Configuration_JSON
template <>
clexmonte::Configuration from_json<clexmonte::Configuration>(
    jsonParser const &json) {
  InputParser<clexmonte::Configuration> parser{json};
  std::stringstream ss;
  ss << "Error reading clexmonte::Configuration from JSON input";
  report_and_throw_if_invalid(parser, CASM::log(),
                              std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - The "basis" must be "prim"
/// - This does not check the validity of the DoF values dimensions
///
/// \ingroup clexmonte_Configuration_JSON
void from_json(clexmonte::Configuration &configuration,
               jsonParser const &json) {
  configuration = from_json<clexmonte::Configuration>(json);
}

/// \@}

}  // namespace CASM
