#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/FixedConfigGenerator.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/crystallography/SymInfo.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"
#include "casm/monte/misc/polymorphic_method_json_io.hh"

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
///     - "fixed": FixedConfigGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "fixed": `parse(InputParser<FixedConfigGenerator> &, ...)`
void parse(
    InputParser<config_generator_type> &parser,
    MethodParserMap<config_generator_type> const &config_generator_methods) {
  parse_polymorphic_method(parser, config_generator_methods);
}

/// \brief Construct FixedConfigGenerator from JSON
///
/// If `configuration` is given, it is used for the initial configuraiton
/// of the Monte Carlo supercell. Otherwise,
/// `transformation_matrix_to_supercell` is used to create the Monte Carlo
/// supercell which is filled with the `motif` configuration.
///
/// Expected format:
/// \code
///   "configuration": object, optional
///       Initial configuration to use for the Monte Carlo supercell. If
///       not given, `motif` must be provided.
///
///   "transformation_matrix_to_supercell": array, shape=3x3
///       Supercell, to be filled with `motif`.
///
///   "motif": object, optional
///       Initial Configuration,
///       which will be copied and tiled into the Monte Carlo supercell.
///       If a perfect tiling can be made by applying factor group operations,
///       a note is printed indicating which operation is applied. A warning
///       is printed if there is no perfect tiling and the `motif` is used
///       without reorientation to fill the supercell imperfectly. If
///       `transformation_matrix_to_supercell` is given but no `motif` is
///       provided, the default configuration is used.
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
void parse(InputParser<FixedConfigGenerator> &parser,
           std::shared_ptr<system_type> const &system) {
  std::unique_ptr<clexmonte::Configuration> configuration;
  if (parser.self.contains("configuration") &&
      !parser.self["configuration"].is_null()) {
    std::unique_ptr<config::Configuration> configuration =
        parser.optional<config::Configuration>("configuration",
                                               *system->supercells);
    parser.value = std::make_unique<FixedConfigGenerator>(*configuration);

  } else if (parser.self.contains("transformation_matrix_to_supercell") &&
             !parser.self["transformation_matrix_to_supercell"].is_null()) {
    auto &log = CASM::log();

    Eigen::Matrix3l T;
    parser.require(T, "transformation_matrix_to_supercell");
    if (!parser.valid()) {
      return;
    }
    std::shared_ptr<config::Supercell const> supercell =
        std::make_shared<config::Supercell const>(system->prim, T);

    if (parser.self.contains("motif") && !parser.self["motif"].is_null()) {
      std::unique_ptr<config::Configuration> motif =
          parser.optional<config::Configuration>("motif", *system->supercells);

      // check if motif can tile into `transformation_matrix_to_supercell`
      auto const &superlattice = supercell->superlattice.superlattice();
      auto const &unit_lattice = motif->supercell->superlattice.superlattice();
      auto const &fg_elements = system->prim->sym_info.factor_group->element;
      double tol = std::max(superlattice.tol(), unit_lattice.tol());
      auto result = is_equivalent_superlattice(superlattice, unit_lattice,
                                               fg_elements.begin(),
                                               fg_elements.end(), tol);
      bool is_equivalent = (result.first != fg_elements.end());
      Index prim_factor_group_index = -1;
      if (is_equivalent) {
        prim_factor_group_index =
            std::distance(fg_elements.begin(), result.first);
        if (prim_factor_group_index != 0) {
          log << "Note: For \"fixed\" configuration generator: " << std::endl;
          log << "Note: `motif` tiles the supercell specified by "
                 "`transformation_matrix_to_supercell` after applying "
                 "prim factor group operation "
              << prim_factor_group_index + 1 << " (indexing from 1)."
              << std::endl;
          log << "Note: Prim factor group operations: (Cartesian)" << std::endl;
          Index i = 1;
          for (auto op : fg_elements) {
            xtal::SymInfo syminfo(op, system->prim->basicstructure->lattice());
            log << "- " << i << ": "
                << to_brief_unicode(syminfo, xtal::SymInfoOptions(CART))
                << std::endl;
            ++i;
          }
          log << std::endl;
        }
      } else {
        log << "Warning: For \"fixed\" configuration generator: " << std::endl;
        log << "Warning: `motif` cannot tile the supercell specified by "
               "`transformation_matrix_to_supercell`. Filling imperfectly.";
        prim_factor_group_index = 0;
      }

      xtal::UnitCell translation(0, 0, 0);
      parser.value = std::make_unique<FixedConfigGenerator>(copy_configuration(
          prim_factor_group_index, translation, *motif, supercell));
    } else {
      log << "Note: For \"fixed\" configuration generator: " << std::endl;
      log << "Note: No \"motif\" parameter. Using default configuration."
          << std::endl;
      parser.value = std::make_unique<FixedConfigGenerator>(
          config::Configuration(supercell));
    }
  } else {
    std::stringstream msg;
    msg << "One of `configuration` or `transformation_matrix_to_supercell` is "
           "required.";
    parser.error.insert(msg.str());
  }
}

}  // namespace clexmonte
}  // namespace CASM
