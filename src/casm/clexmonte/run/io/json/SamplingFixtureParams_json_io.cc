#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"

#include "casm/monte/MethodLog.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"
#include "casm/monte/run_management/SamplingFixture.hh"
#include "casm/monte/run_management/io/ResultsIO.hh"
#include "casm/monte/run_management/io/json/ResultsIO_json_io.hh"
#include "casm/monte/run_management/io/json/jsonResultsIO_impl.hh"
#include "casm/monte/sampling/io/json/SamplingParams_json_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct sampling_fixture_params_type from JSON
///
///
/// \code
/// {
///     "sampling": <monte::SamplingParams>
///         Options controlling which quantities are sampled and how often
///         sampling is performed.
///     "completion_check": <monte::CompletionCheck>
///         Controls when a sampling fixture is complete. Options include
///         convergence of sampled quantiies, min/max number of samples, min/
///         max number of passes, etc.
///     "analysis":
///         "functions": array of str (default=[])
///             Names of analysis functions to use to evaluate results with.
///             Standard options include the following (others may be included):
///
///             - "heat_capacity": Heat capacity
///             - "mol_susc": Chemical susceptibility (mol_composition)
///             - "param_susc": Chemical susceptibility (param_composition)
///             - "mol_thermocalc_susc": Thermo-chemical susceptibility
///               (mol_composition)
///             - "param_thermocalc_susc": Thermo-Chemical susceptibility
///               (param_composition)
///
///             Unless otherwise noted, assume per unitcell properties.
///
///     "results_io": <monte::ResultsIO> = null
///         Options controlling sampling fixture results output.
///     "log": (optional)
///         "file": str (default="status.json")
///             Provide the path where a log file should be written.
///         "frequency_in_s": number (default=600.0)
///             How often the log file should be written, in seconds.
/// }
/// \endcode
void parse(InputParser<sampling_fixture_params_type> &parser, std::string label,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions,
           MethodParserMap<results_io_type> const &results_io_methods,
           bool time_sampling_allowed) {
  // Read sampling params
  std::set<std::string> sampling_function_names;
  for (auto const &element : sampling_functions) {
    sampling_function_names.insert(element.first);
  }
  auto sampling_params_subparser = parser.subparse<monte::SamplingParams>(
      "sampling", sampling_function_names, time_sampling_allowed);
  if (!parser.valid()) {
    return;
  }
  monte::SamplingParams const &sampling_params =
      *sampling_params_subparser->value;
  std::map<std::string, state_sampling_function_type>
      selected_sampling_functions;
  for (auto const &name : sampling_params.sampler_names) {
    selected_sampling_functions.emplace(name, sampling_functions.at(name));
  }

  // Read completion check params
  auto completion_check_params_subparser =
      parser.subparse<monte::CompletionCheckParams<statistics_type>>(
          "completion_check", sampling_functions);

  // Read analysis functions
  std::vector<std::string> function_names;
  fs::path functions_path = fs::path("analysis") / "functions";
  parser.optional(function_names, functions_path);

  std::map<std::string, results_analysis_function_type>
      selected_analysis_functions;
  for (auto const &name : function_names) {
    auto it = analysis_functions.find(name);
    if (it != analysis_functions.end()) {
      selected_analysis_functions.insert(*it);
    } else {
      std::stringstream msg;
      msg << "Error: function '" << name << "' not recognized";
      parser.insert_error(functions_path, msg.str());
    }
  }

  // Construct results I/O instance
  auto results_io_subparser =
      parser.subparse_if<results_io_type>("results_io", results_io_methods);

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
    parser.value = std::make_unique<sampling_fixture_params_type>(
        label, selected_sampling_functions, selected_analysis_functions,
        sampling_params, *completion_check_params_subparser->value,
        std::move(results_io_subparser->value), method_log);
  }
}

}  // namespace clexmonte
}  // namespace CASM
