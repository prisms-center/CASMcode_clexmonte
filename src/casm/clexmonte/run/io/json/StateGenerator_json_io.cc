#include "casm/clexmonte/run/io/json/StateGenerator_json_io_impl.hh"
#include "casm/monte/misc/polymorphic_method_json_io.hh"

namespace CASM {
namespace clexmonte {

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
/// Expected JSON:
///   method: string (required)
///     The name of the chosen state generation method. Currently, the only
///     option is:
///     - "incremental": IncrementalConditionsStateGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "incremental":
///           `parse(InputParser<incremental_state_generator_type> &, ...)`
///
void parse(
    InputParser<state_generator_type> &parser,
    MethodParserMap<state_generator_type> const &state_generator_methods) {
  parse_polymorphic_method(parser, state_generator_methods);
}

}  // namespace clexmonte
}  // namespace CASM
