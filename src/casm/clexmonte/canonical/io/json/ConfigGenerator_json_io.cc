#include "casm/clexmonte/canonical/io/json/ConfigGenerator_json_io.hh"

#include "casm/clexmonte/clex/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct ConfigGenerator from JSON
void parse(InputParser<config_generator_type> &parser,
           std::shared_ptr<system_type> const &system_data) {
  PolymorphicParserFactory<config_generator_type> f;
  parse_polymorphic_method(
      parser, {f.make<fixed_config_generator_type>("fixed", system_data)});
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
