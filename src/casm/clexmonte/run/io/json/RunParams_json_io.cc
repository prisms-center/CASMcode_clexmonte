#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/run/io/json/ResultsIO_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/sampling/io/json/SamplingParams_json_io.hh"
//#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse canonical Monte Carlo input file
///
/// Input file summary:
/// \code
/// {
///     "system": <clexmonte::System>
///         Species the path to a "system" input file.
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
/// }
/// \endcode
void parse(
    InputParser<RunParams> &parser, std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions) {
  // Construct state generator
  auto state_generator_subparser = parser.subparse<state_generator_type>(
      "state_generation", system, sampling_functions);

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
  auto results_io_subparser = parser.subparse<results_io_type>(
      "results_io", sampling_functions, analysis_functions);

  if (parser.valid()) {
    parser.value =
        std::make_unique<RunParams>(sampling_functions, analysis_functions,
                                    std::move(state_generator_subparser->value),
                                    *sampling_params_subparser->value,
                                    *completion_check_params_subparser->value,
                                    std::move(results_io_subparser->value));
  }
}

}  // namespace clexmonte
}  // namespace CASM
