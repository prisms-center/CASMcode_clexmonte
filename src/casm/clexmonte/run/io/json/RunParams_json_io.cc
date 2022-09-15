#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/run/io/json/ResultsIO_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/results/io/json/jsonResultsIO_impl.hh"
#include "casm/monte/sampling/io/json/SamplingParams_json_io.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

MethodParserMap<config_generator_type> standard_config_generator_methods(
    std::shared_ptr<system_type> const &system) {
  MethodParserFactory<config_generator_type> cf;
  MethodParserMap<config_generator_type> config_generator_methods;
  config_generator_methods.insert(
      cf.make<monte::FixedConfigGenerator<config_type>>("fixed", system));
  return config_generator_methods;
}

MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods) {
  MethodParserFactory<state_generator_type> sf;
  MethodParserMap<state_generator_type> state_generator_methods;
  state_generator_methods.insert(
      sf.make<monte::IncrementalConditionsStateGenerator<Configuration>>(
          "incremental", system, sampling_functions, config_generator_methods));
  return state_generator_methods;
}

MethodParserMap<results_io_type> standard_results_io_methods(
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions) {
  MethodParserFactory<results_io_type> f;
  MethodParserMap<results_io_type> results_io_methods;
  results_io_methods.insert(f.template make<monte::jsonResultsIO<config_type>>(
      "json", sampling_functions, analysis_functions));
  return results_io_methods;
}

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
///     "log": (optional)
///       "file": str (default="status.json")
///         Provide the path where a log file should be written.
///       "frequency_in_s": number (default=600.0)
///         How often the log file should be written, in seconds.
/// }
/// \endcode
void parse(
    InputParser<RunParams> &parser, std::shared_ptr<system_type> const &system,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions,
    MethodParserMap<state_generator_type> const &state_generator_methods,
    MethodParserMap<results_io_type> const &results_io_methods) {
  // Construct state generator
  auto state_generator_subparser = parser.subparse<state_generator_type>(
      "state_generation", state_generator_methods);

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
      parser.subparse<results_io_type>("results_io", results_io_methods);

  // Method log
  monte::MethodLog method_log;
  if (parser.self.contains("log")) {
    std::string log_file = "status.json";
    parser.optional(log_file, fs::path("log") / "file");
    double log_frequency = 600.0;
    parser.optional(log_frequency, fs::path("log") / "frequency_in_s");

    method_log.log_frequency = log_frequency;
    method_log.logfile_path = fs::path(log_file);
    method_log.reset();
  }

  if (parser.valid()) {
    parser.value = std::make_unique<RunParams>(
        sampling_functions, analysis_functions,
        std::move(state_generator_subparser->value),
        *sampling_params_subparser->value,
        *completion_check_params_subparser->value,
        std::move(results_io_subparser->value), method_log);
  }
}

}  // namespace clexmonte
}  // namespace CASM
