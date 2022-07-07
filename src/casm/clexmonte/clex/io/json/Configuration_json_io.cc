#include "casm/clexmonte/clex/io/json/Configuration_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/Configuration.hh"
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
///
/// ### Examples
///
/// #### Example 1) Configuration with occupation, displacement, and strain DoF
/// ```
/// {
///   "transformation_matrix_to_super": [],
///   "dof": {}
/// }
/// ```

/// \brief Write clexmonte::Configuration to JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values are written in the
///   basis in which they are provided
///
/// \ingroup clexmonte_Configuration_JSON
jsonParser &to_json(clexmonte::Configuration const &configuration,
                    jsonParser &json) {
  json["transformation_matrix_to_supercell"] =
      get_transformation_matrix_to_super(configuration);
  json["dof"] = get_dof_values(configuration);
  return json;
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
///
/// \ingroup clexmonte_Configuration_JSON
template <>
clexmonte::Configuration from_json<clexmonte::Configuration>(
    jsonParser const &json) {
  ParentInputParser parser{json};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");

  clexulator::ConfigDoFValues dof_values;
  parser.require(dof_values, "dof");

  auto &log = CASM::log();
  std::runtime_error error_if_invalid{
      "Error reading clexmonte::Configuration from JSON input"};
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  return clexmonte::Configuration(T, dof_values);
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
///
/// \ingroup clexmonte_Configuration_JSON
void from_json(clexmonte::Configuration &configuration,
               jsonParser const &json) {
  configuration = from_json<clexmonte::Configuration>(json);
}

/// \@}

}  // namespace CASM
