#include "casm/clexmonte/canonical/io/json/StateGenerator_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/io/json/Conditions_json_io.hh"
#include "casm/clexmonte/canonical/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/clex/Configuration.hh"
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
///     - "incremental"
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "incremental": `parse_canonical_state_generator`
///
void parse(
    InputParser<state_generator_type> &parser,
    std::shared_ptr<system_type> const &system_data,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    canonical_tag tag) {
  PolymorphicParserFactory<state_generator_type> f;
  parse_polymorphic_method(
      parser, {f.make<incremental_state_generator_type>(
                  "incremental", system_data, sampling_functions, tag)});
}

/// \brief Construct IncrementalConditionsStateGenerator from JSON
///
/// The "incremental" state generation method generates a series of N states,
/// starting at a particular set of conditions and then incrementing by a fixed
/// amount the conditions values for each subsequent step. It allows users to
/// specify:
/// - independent conditions, which are explicitly specified by a choice
///   of initial value and incremental value;
/// - a `ConfigGenerator` method, which generates a configuration as a function
///   of the independent conditions, the final states of previous runs, and the
///   calculation results of previous runs;
/// - and dependent conditions, which are specified as a function of the
///   independent conditions and generated configuration.
///
/// For canonical Monte Carlo calculations, a common use case is to specify
/// the temperature range independently via `initial_conditions` and
/// `conditions_increment`, while letting the composition be a dependent
/// condition of the generated configuration. This results in calculations at
/// fixed composition and varying temperature, which is often of interest, and
/// relieves the user of having to calculate the composition of the initial
/// configuration and set it manually.
///
/// Expected:
///   initial_configuration: ConfigGenerator
///     Specifies how to generate an initial configuration given the
///     current conditions and previous runs.
///
///     method: string (required)
///       Choice of ConfigGenerator method. Currently, the only option is:
///       - "fixed"
///
///     kwargs: object (optional)
///        Options for chosen ConfigGenerator method.
///
///        For the "fixed" method, this object has the same format as a
///        Configuration. The attributes are:
///
///          transformation_matrix_to_super: 3x3 matrix of int
///            Matrix, T, specifiying the supercell lattice vectors in terms of
///            the primitive cell vectors, according to `S = P * T`, where the
///            columns of P are the primitive lattice vectors (a,b,c) and the
///            columns of S are the supercell lattice vectors.
///
///          dof: ConfigDoFValues object (default=zeros configuration)
///            Degree of freedom (DoF) values in the supercell defined by
///            `transformation_matrix_to_super`, and expressed in the *standard
///            basis*. See the ConfigDoFValues format for details.
///
///            If the value is `null` or not provided, the default
///            is that all DoF values are set equal to zero. For occupation
///            DoF, the default value of zero means the first type listed on
///            the corresponding sublattice in the `"prim"`.
///
///   initial_conditions: object
///     Conditions for the initial state. For canonical Monte Carlo
///     calculations, "temperature" is required and composition (using "comp_n"
///     or "comp_x") must be specified for "initial_conditions" or listed in
///     "dependent_conditions". May include:
///
///       "temperature": number (required)
///         Temperature in K.
///
///       "comp_n": array of number or dict (optional)
///         Composition in number per primitive cell. May be:
///
///         - An array of number, specifying the number of each component per
///           primitive cell, interpreted using the order of the `"components"`
///           specified for composition axes. The size must match the number of
///           components.
///
///         - A dict, where the keys are the component names, and values are
///           the number of that component per primitive cell. All components in
///           the system must be included.
///
///       "comp_x": array of number or dict (optional)
///         Parametric composition, in terms of the chosen composition axes. May
///         be:
///
///         - An array of number, corresponding to
///           `[comp_x(a), comp_x(b), ...]`. The size must match the number of
///           composition axes.
///
///         - A dict, where the keys are the axes names ("a", "b", etc.), and
///           values are the corresponding parametric composition value.
///           All composition axes must be included.
///
///   conditions_increment: object (required)
///     Amount to increment the independent conditions for each subsequent
///     state. All conditions listed for `"initial_conditions"` must be
///     specified here, even if the increment is zero valued or only a single
///     state will be generated.
///
///   n_state: integer (required)
///     Total number of states to generate. Includes the inital state.
///
///   dependent_runs: bool (optional, default=true)
///     If true, only use the ConfigGenerator specified by
///     "initial_configuration" for the first state. For subsequent states at
///     new conditions, use the configuration of the final state for the
///     initial configuration at the next condistions. Choosing `true` tends to
///     result in smoother calculation results from condition to condition and
///     more hysteresis.
///
///     If false, then always use the ConfigGenerator to determine the initial
///     configuration. Choosing `false` tends to result in noisier calculation
///     results from condition to condition and less hysteresis.
///
///   dependent_conditions: Array of string (optional, default=[])
///     Names of sampling functions specifying conditions that should be
///     fixed at the value the results from evaluating the state generated
///     by the choice of initial configuration and increment conditions.
///
///     For canonical Monte Carlo calculations, a common use case is to specify
///     the temperature range independently via `initial_conditions` and
///     `conditions_increment` while fixing the composition of the initial
///     configuration by using `"dependent_conditions": ["comp_n"]`. This
///     relieves the user of having to calculate the composition of the initial
///     configuration and set it in `initial_conditions` manually.
///
void parse(
    InputParser<incremental_state_generator_type> &parser,
    std::shared_ptr<system_type> const &system_data,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions) {
  /// Helps select otherwise overloaded `parse` methods
  canonical_tag tag;

  /// Parse "initial_configuration"
  auto config_generator_subparser = parser.subparse<config_generator_type>(
      "initial_configuration", system_data, tag);

  /// Parse "initial_conditions"
  auto initial_conditions_subparser =
      parser.subparse_with<monte::VectorValueMap>(
          canonical::parse_conditions, "initial_conditions",
          get_composition_converter(*system_data), tag);

  /// Parse "conditions_increment"
  auto conditions_increment_subparser =
      parser.subparse_with<monte::VectorValueMap>(
          canonical::parse_conditions_increment, "conditions_increment",
          get_composition_converter(*system_data), tag);

  /// Parse "dependent_conditions"
  std::vector<std::string> dependent_conditions_names;
  parser.optional(dependent_conditions_names, "dependent_conditions");
  monte::StateSamplingFunctionMap<config_type> dependent_conditions;
  for (auto const &name : dependent_conditions_names) {
    auto it = sampling_functions.find(name);
    if (it == sampling_functions.end()) {
      std::stringstream msg;
      msg << "Error in \"dependent_conditions\": Not a valid sampling function "
             "name: \""
          << name << "\"";
      parser.insert_error("dependent_conditions", msg.str());
      continue;
    }
    dependent_conditions.emplace(*it);
  }

  /// Parse "n_states"
  Index n_states;
  parser.require(n_states, "n_states");

  /// Parse "dependent_runs"
  bool dependent_runs = true;
  parser.optional(dependent_runs, "dependent_runs");

  if (parser.valid()) {
    parser.value = std::make_unique<incremental_state_generator_type>(
        std::move(config_generator_subparser->value),
        *initial_conditions_subparser->value,
        *conditions_increment_subparser->value, n_states, dependent_runs,
        dependent_conditions);
  }
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
