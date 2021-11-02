#include "casm/clexmonte/canonical/io/json/ConfigGenerator_json_io.hh"

#include "casm/clexmonte/clex/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

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
void parse(InputParser<config_generator_type> &parser,
           std::shared_ptr<system_type> const &system_data, canonical_tag tag) {
  PolymorphicParserFactory<config_generator_type> f;
  parse_polymorphic_method(
      parser, {f.make<fixed_config_generator_type>("fixed", system_data)});
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
