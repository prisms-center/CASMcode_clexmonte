#include "casm/clexmonte/canonical/io/json/InputData_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/io/InputData.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/results/io/json/ResultsIO_json_io_impl.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/state/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/io/json/parse_conditions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/sampling/io/json/SamplingParams_json_io.hh"
#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct conditions (monte::ValueMap) from JSON
///
/// Parses canonical Monte Carlo conditions from JSON. If successfully parsed,
/// `parser->value` will contain a monte::ValueMap with:
/// - scalar_values["temperature"]
/// - vector_values["mol_composition"]: (size = system components size)
///
/// Expected:
///
///   "temperature": number (required)
///     Temperature in K.
///
///   "mol_composition": dict (optional)
///     Composition in number per primitive cell. A dict, where the keys are
///     the component names, and values are the number of that component per
///     primitive cell. All components in the system must be included.
///
///   "param_composition": array of number or dict (optional)
///     Parametric composition, in terms of the chosen composition axes. Will
///     be converted to `"mol_composition"`. A dict, where the keys are the axes
///     names
///     ("a", "b", etc.), and values are the corresponding parametric
///     composition value. All composition axes must be included.
///
///
void parse_conditions(InputParser<monte::ValueMap> &parser,
                      std::shared_ptr<system_type> const &system,
                      canonical_tag tag) {
  parser.value = std::make_unique<monte::ValueMap>();
  parse_temperature(parser);
  parse_mol_composition(parser, system);
}

/// \brief Construct conditions increment (monte::ValueMap) from JSON
///
/// Parses canonical Monte Carlo conditions increments from JSON. If
/// successfully parsed, `parser->value` will contain a monte::ValueMap
/// with:
/// - "temperature": (size 1)
/// - "mol_composition": (size = system components size)
///
/// The expected JSON format is the same as documented for `parse_conditions`,
/// but values are interpreted as increments.
void parse_conditions_increment(InputParser<monte::ValueMap> &parser,
                                std::shared_ptr<system_type> const &system,
                                canonical_tag tag) {
  parser.value = std::make_unique<monte::ValueMap>();
  parse_temperature(parser);
  parse_mol_composition_increment(parser, system);
}

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
           std::shared_ptr<system_type> const &system, canonical_tag tag) {
  PolymorphicParserFactory<config_generator_type> f;
  parse_polymorphic_method(
      parser, {f.make<fixed_config_generator_type>("fixed", system)});
}

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
///     - "incremental": monte::IncrementalConditionsStateGenerator
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "incremental":
///           `parse(InputParser<incremental_state_generator_type> &, ...)`
///
void parse(
    InputParser<state_generator_type> &parser,
    std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    canonical_tag tag) {
  PolymorphicParserFactory<state_generator_type> f;
  parse_polymorphic_method(parser,
                           {f.make<incremental_state_generator_type>(
                               "incremental", system, sampling_functions,
                               canonical::parse_conditions,
                               canonical::parse_conditions_increment, tag)});
}

/// \brief Parse canonical Monte Carlo input file
///
/// Input file summary:
/// \code
/// {
///   "method": "canonical"
///   "kwargs": {
///     "system": {
///       "prim": <xtal::BasicStructure or file path>
///           Specifies the primitive crystal structure and allowed DoF. Must
///           be the prim used to generate the cluster expansion.
///       "composition_axes": <composition::CompositionConverter>
///           Specifies composition axes
///       "formation_energy": <Clex>
///           Input specifies a CASM cluster expansion basis set source file,
///           coefficients, and compilation settings.
///     },
///     "state_generation": <monte::StateGenerator>
///         Specifies a "path" of input states at which to run Monte Carlo
///         calculations. Each state is an initial configuration and set of
///         thermodynamic conditions (temperature, chemical potential,
///         composition, etc.).
///     "random_number_generator": <monte::RandomNumberGenerator>
///         (Future) Options controlling the random number generator.
///     "sampling": <monte::SamplingParams>
///         Options controlling which quantities are sampled and how often
///         sampling is performed.
///     "completion_check": <monte::CompletionCheck>
///         Controls when a single Monte Carlo run is complete. Options include
///         convergence of sampled quantiies, min/max number of samples, min/
///         max number of passes, etc.
///     "results_io": <monte::ResultsIO>
///         Options controlling results output.
///   }
/// }
void parse(InputParser<InputData> &parser) {
  // Parse canonical MC calculation data. Includes input:
  // - "prim"
  // - "composition_axes"
  // - "formation_energy"
  auto system_subparser = parser.subparse<system_type>("system");
  if (!system_subparser->valid()) {
    return;
  }
  std::shared_ptr<system_type> system = std::move(system_subparser->value);

  // Make state sampling functions, with current supercell-specific info
  monte::StateSamplingFunctionMap<config_type> sampling_functions =
      make_sampling_functions(system, canonical_tag());

  // Construct state generator
  auto state_generator_subparser = parser.subparse<state_generator_type>(
      "state_generation", system, sampling_functions, canonical_tag());

  // Read sampling params
  std::set<std::string> sampling_function_names;
  for (auto const &element : sampling_functions) {
    sampling_function_names.insert(element.first);
  }
  bool time_sampling_allowed = false;
  auto sampling_params_subparser = parser.subparse<monte::SamplingParams>(
      "sampling", sampling_function_names, time_sampling_allowed);

  // Read completion check params
  auto completion_check_params_subparser =
      parser.subparse<monte::CompletionCheckParams>("completion_check",
                                                    sampling_functions);

  // Construct results I/O instance
  auto results_io_subparser =
      parser.subparse<results_io_type>("results_io", sampling_functions);

  // Construct random number generator
  MTRand random_number_generator;

  if (parser.valid()) {
    parser.value = std::make_unique<InputData>(
        system, std::move(state_generator_subparser->value), sampling_functions,
        *sampling_params_subparser->value,
        *completion_check_params_subparser->value,
        std::move(results_io_subparser->value), random_number_generator);
  }
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
