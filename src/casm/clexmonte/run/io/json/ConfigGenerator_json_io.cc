#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/run_management/FixedConfigGenerator.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct ConfigGenerator from JSON
///
/// A configuration generation method generates a configuration given a set of
/// conditions and results from previous runs. It may be a way to customize a
/// state generation method.
///
/// Expected:
///   method: string (required)
///     The name of the chosen config generation method. Currently, the only
///     option is:
///     - "fixed": monte::FixedConfigGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "fixed": `parse(InputParser<monte::FixedConfigGenerator> &, ...)`
void parse(
    InputParser<config_generator_type> &parser,
    MethodParserMap<config_generator_type> const &config_generator_methods) {
  parse_polymorphic_method(parser, config_generator_methods);
}

/// \brief Construct FixedConfigGenerator from JSON
///
/// Expected format:
/// \code
///   "transformation_matrix_to_supercell": array, shape=3x3
///       Supercell
///
///   "dof": object, optional
///       Initial ConfigDoFValues, in standard basis, for the Monte
///       Carlo supercell. If no initial given, the default
///       configuration is used.
///
///   "motif": object, optional
///       Initial Configuration, with ConfigDoFValues in standard basis,
///       which will be copied and tiled into the Monte Carlo supercell.
///       There is no warning if the tiling is not perfect. If no
///       initial given, the default configuration is used.
///
/// \endcode
///
///
/// Requires:
/// - `Configuration from_standard_values(
///        system_type const &system,
///        Configuration const &configuration)`
/// - `Configuration make_default_configuration(
///        system_type const &system,
///        Eigen::Matrix3l const &transformation_matrix_to_super)`
void parse(InputParser<monte::FixedConfigGenerator<Configuration>> &parser,
           std::shared_ptr<system_type> const &system) {
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  if (!parser.valid()) {
    return;
  }

  std::unique_ptr<clexmonte::Configuration> configuration;
  if (parser.self.contains("dof") && !parser.self["dof"].is_null()) {
    // TODO: validation of dof types and dimensions?
    // Note: expect "dof" to be in standard basis
    std::unique_ptr<clexulator::ConfigDoFValues> standard_dof_values =
        parser.optional<clexulator::ConfigDoFValues>("dof");
    configuration = std::make_unique<Configuration>(T, *standard_dof_values);
  } else if (parser.self.contains("motif") && !parser.self["motif"].is_null()) {
    // Note: expect motif dof to be in standard basis
    std::unique_ptr<config::Configuration> motif =
        parser.optional<config::Configuration>("motif", system->prim);
    std::shared_ptr<config::Supercell const> supercell =
        std::make_shared<config::Supercell const>(system->prim, T);
    config::Configuration tmp = copy_configuration(*motif, supercell);
    configuration = std::make_unique<Configuration>(T, tmp.dof_values);
  }

  if (parser.valid()) {
    if (configuration != nullptr) {
      clexmonte::Configuration tmp =
          from_standard_values(*system, *configuration);
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(tmp);
    } else {
      parser.value =
          notstd::make_unique<monte::FixedConfigGenerator<Configuration>>(
              make_default_configuration(*system, T));
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM
