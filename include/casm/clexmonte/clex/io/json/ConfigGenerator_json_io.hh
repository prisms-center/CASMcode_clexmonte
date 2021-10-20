#ifndef CASM_clexmonte_clex_ConfigGenerator_json_io
#define CASM_clexmonte_clex_ConfigGenerator_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct monte::FixedConfigGenerator from JSON
template <typename SystemType>
void parse(InputParser<monte::FixedConfigGenerator<Configuration>> &parser,
           std::shared_ptr<SystemType> const &system_data);

// --- Inline implementations ---

/// \brief Construct FixedConfigGenerator from JSON
///
/// Requires:
/// - `Configuration from_standard_values(
///        SystemType const &system_data,
///        Configuration const &configuration)`
/// - `Configuration make_default_configuration(
///        SystemType const &system_data,
///        Eigen::Matrix3l const &transformation_matrix_to_super)`
template <typename SystemType>
void parse(InputParser<monte::FixedConfigGenerator<Configuration>> &parser,
           std::shared_ptr<SystemType> const &system_data) {
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_super");
  if (!parser.valid()) {
    return;
  }

  // TODO: validation of dof types and dimensions?
  // Note: expect "dof" to be in standard basis
  std::unique_ptr<clexulator::ConfigDoFValues> standard_dof_values =
      parser.optional<clexulator::ConfigDoFValues>("dof");

  if (parser.valid()) {
    if (standard_dof_values != nullptr) {
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(
              from_standard_values(*system_data,
                                   Configuration(T, *standard_dof_values)));
    } else {
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(
              make_default_configuration(*system_data, T));
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
