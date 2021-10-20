#include "casm/clexmonte/io/json/Configuration_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/Configuration.hh"

namespace CASM {

/// \brief Write clexmonte::Configuration to JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values are written in the
///   basis in which they are provided
jsonParser &to_json(clexmonte::Configuration const &configuration,
                    jsonParser &json) {
  json["transformation_matrix_to_supercell"] =
      configuration.transformation_matrix_to_super();
  json["dof"] = configuration.dof_values;
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
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

  return Configuration(T, dof_values);
}

/// \brief Read clexmonte::Configuration from JSON
///
/// Notes:
/// - This does not convert the DoF values basis, values stay in the basis in
/// which they are provided
/// - This does not check the validity of the DoF values dimensions
void from_json(clexmonte::Configuration &configuration,
               jsonParser const &json) {
  f = from_json<clexmonte::Configuration>(json);
}

}  // namespace CASM
