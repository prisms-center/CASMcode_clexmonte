#ifndef CASM_clexmonte_run_StateGenerator_json_io_impl
#define CASM_clexmonte_run_StateGenerator_json_io_impl

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/IncrementalConditionsStateGenerator.hh"
#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/parse_conditions.hh"
// #include "casm/monte/misc/polymorphic_method_json_io.hh"

namespace CASM {
namespace clexmonte {

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
///          configuration: object, optional
///              Initial configuration to use for the Monte Carlo supercell. If
///              not given, `motif` must be provided.
///
///          transformation_matrix_to_supercell: array, shape=3x3
///              Supercell, to be filled with `motif`.
///
///          motif: object, optional
///              Initial Configuration, which will be copied and tiled into
///              the Monte Carlo supercell. If a perfect tiling can be made
///              by applying factor group operations, a note is printed
///              indicating which operation is applied. A warning is printed
///              if there is no perfect tiling and the `motif` is used
///              without reorientation to fill the supercell imperfectly. If
///              `transformation_matrix_to_supercell` is given but no `motif` is
///              provided, the default configuration is used.
///
///   initial_conditions: object
///     Conditions for the initial state. For canonical Monte Carlo
///     calculations, "temperature" is required and composition (using
///     "mol_composition" or "param_composition") must be specified for
///     "initial_conditions" or set via "modifiers". May include:
///
///       "temperature": number (required)
///         Temperature in K.
///
///       "mol_composition": array of number or dict (optional)
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
///       "param_composition": array of number or dict (optional)
///         Parametric composition, in terms of the chosen composition axes. May
///         be:
///
///         - An array of number, specifying the parametric composition along
///           each axis (i.e. `[a, b, ...]`). The size must match the number
///           of composition axes.
///
///         - A dict, where the keys are the axes names ("a", "b", etc.), and
///           values are the corresponding parametric composition value.
///           All composition axes must be included.
///
///       "param_chem_pot": array of number or dict (optional)
///         Parametric chemical potential, i.e. the chemical potential conjugate
///         to the chosen composition axes. May be:
///
///         - An array of number, specifying the parametric chemical potential
///           for each axis (i.e. `[a, b, ...]`). The size must match the number
///           of composition axes.
///
///         - A dict, where the keys are the axes names ("a", "b", etc.), and
///           values are the corresponding parametric chemical potential
///           values. All composition axes must be included.
///
///
///   conditions_increment: object (required)
///     Amount to increment the independent conditions for each subsequent
///     state. All conditions listed for `"initial_conditions"` must be
///     specified here, even if the increment is zero valued or only a single
///     state will be generated.
///
///   n_states: integer (required)
///     Total number of states to generate. Includes the initial state.
///
///   dependent_runs: bool (optional, default=true)
///     If true, only use the ConfigGenerator specified by
///     "initial_configuration" for the first state. For subsequent states at
///     new conditions, use the configuration of the final state for the
///     initial configuration at the next conditions. Choosing `true` tends to
///     result in smoother calculation results from condition to condition and
///     more hysteresis.
///
///     If false, then always use the ConfigGenerator to determine the initial
///     configuration. Choosing `false` tends to result in noisier calculation
///     results from condition to condition and less hysteresis.
///
///   completed_runs: dict (optional)
///     Controls storage and output of completed runs data used for handling
///     restarts. One of "save_all_final_states" or "save_last_final_state"
///     must be true for the "dependent_runs" option. May include:
///
///       "save_all_initial_states": bool = false
///         If true, save initial states for analysis, output, or state
///         generation.
///
///       "save_all_final_states": bool = false
///         If true, save final states for analysis, output, or state
///         generation.
///
///       "save_last_final_state": bool = true
///         If true, save final state for last run for analysis, output
///         or state generation.
///
///       "write_initial_states": bool = false
///         If true, write saved initial states to completed_runs.json.
///
///       "write_final_states": bool = false
///         If true, write saved final states to completed_runs.json.
///
///       "output_dir": str = ""
///         If not empty, name of a directory in which to write
///         completed_runs.json. If empty, completed_runs.json is not
///         written, which means restarts are not possible.
///
///   modifiers: Array of string (optional, default=[])
///     Names of functions that should be used to modify the state generated
///     by the choice of initial configuration and conditions increment. The
///     modifiers may in general modify either the configuration or conditions
///     of the state.
///
///     For canonical Monte Carlo calculations, a common use case is to specify
///     the temperature range independently via `initial_conditions` and
///     `conditions_increment` while fixing the composition of the initial
///     configuration by using `"modifiers": ["set_mol_composition"]`.
///     This relieves the user of having to calculate the composition of the
///     initial configuration and set it in `initial_conditions` manually.
///
template <typename ConditionsType>
void parse(InputParser<IncrementalConditionsStateGenerator> &parser,
           std::shared_ptr<system_type> const &system,
           StateModifyingFunctionMap const &modifying_functions,
           MethodParserMap<config_generator_type> config_generator_methods,
           ConditionsType const *ptr) {
  /// Parse "initial_configuration"
  auto config_generator_subparser = parser.subparse<config_generator_type>(
      "initial_configuration", config_generator_methods);

  /// Parse "initial_conditions"
  bool is_increment = false;
  auto initial_conditions_subparser = parser.subparse<ConditionsType>(
      "initial_conditions", system, is_increment);

  /// Parse "conditions_increment"
  is_increment = true;
  //  auto conditions_increment_subparser =
  //  parser.subparse_with<monte::ValueMap>(
  //      parse_conditions, "conditions_increment", system, is_increment);
  auto conditions_increment_subparser = parser.subparse<ConditionsType>(
      "conditions_increment", system, is_increment);

  /// Parse "modifiers"
  std::vector<std::string> modifier_names;
  parser.optional(modifier_names, "modifiers");
  std::vector<StateModifyingFunction> selected_modifiers;
  for (auto const &name : modifier_names) {
    auto it = modifying_functions.find(name);
    if (it == modifying_functions.end()) {
      std::stringstream msg;
      msg << "Error in \"modifiers\": Not a valid function "
             "name: \""
          << name << "\"";
      parser.insert_error("modifiers", msg.str());
      continue;
    }
    selected_modifiers.push_back(it->second);
  }

  /// Parse "n_states"
  Index n_states;
  parser.require(n_states, "n_states");

  /// Parse "dependent_runs"
  bool dependent_runs = true;
  parser.optional(dependent_runs, "dependent_runs");

  /// Parse "completed_runs"
  RunDataOutputParams output_params;
  parser.optional(output_params, "completed_runs");

  if (parser.valid()) {
    parser.value = std::make_unique<IncrementalConditionsStateGenerator>(
        system, output_params, std::move(config_generator_subparser->value),
        initial_conditions_subparser->value->to_value_map(false),
        conditions_increment_subparser->value->to_value_map(true), n_states,
        dependent_runs, selected_modifiers);
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
