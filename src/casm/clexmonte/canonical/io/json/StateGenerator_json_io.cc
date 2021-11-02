#include "casm/clexmonte/canonical/io/json/StateGenerator_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/io/json/Conditions_json_io.hh"
#include "casm/clexmonte/canonical/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/clex/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct StateGenerator from JSON
///
/// A state generation method generates the initial state for each run in a
/// series of Monte Carlo calculation. A state consists of:
/// - a configuration, the choice of periodic supercell lattice vectors and the
/// values of degrees of freedom (DoF) in that supercell along with any global
/// DoF.
/// - a set of thermodynamic conditions, which control the statistical ensemble
/// used. In general, this may include quantities such as temperature, chemical
/// potential, composition, pressure, volume, strain, magnetic field, etc.
/// depending on the type of calculation.
///
/// Expected:
///   method: string (required)
///     The name of the chosen state generation method. Currently, the only
///     option is:
///     - "incremental": monte::IncrementalConditionsStateGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "incremental":
///           `parse(InputParser<incremental_state_generator_type> &, ...)`
///
void parse(
    InputParser<state_generator_type> &parser,
    std::shared_ptr<system_type> const &system_data,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    canonical_tag tag) {
  PolymorphicParserFactory<state_generator_type> f;
  parse_polymorphic_method(parser,
                           {f.make<incremental_state_generator_type>(
                               "incremental", system_data, sampling_functions,
                               canonical::parse_conditions,
                               canonical::parse_conditions_increment, tag)});
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
