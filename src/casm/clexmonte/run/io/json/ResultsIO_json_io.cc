#include "casm/clexmonte/run/io/json/ResultsIO_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/monte/results/io/json/jsonResultsIO_impl.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct ResultsIO from JSON
///
/// A ResultsIO method implements reading and writing Monte Carlo output.
///
/// Expected:
///   method: string (required)
///     The name of the chosen state generation method. Currently, the only
///     option is:
///     - "json": monte::jsonResultsIO
///
///   kwargs: dict (optional, default={})
///     Method-specific options. See documentation for particular methods:
///     - "json": `parse(monte::jsonResultsIO<config_type> &)`
///
void parse(InputParser<results_io_type> &parser,
           MethodParserMap<results_io_type> const &results_io_methods) {
  parse_polymorphic_method(parser, results_io_methods);
}

/// \brief Construct jsonResultsIO from JSON
///
/// The "json" results IO method writes results to JSON files:
/// - summary.json: Summarizes results from each inidividual run in arrays
/// - run.<index>/trajectory.json: Optional output file (one per each
///   individual run) is a JSON array of the configuration at the time a sample
///   was taken.
/// - run.<index>/observations.json: Optional output file (one per each
///   individual run) contains all sampled data.
///
/// Expected:
///   output_dir: string (required)
///     Specifies the directory where results should be written.
///
///   write_observations: bool (default=false)
///     If true, write an `"observations.json"` file for each individual run.
///
///   write_trajectory: bool (default=false)
///     If true, write an `"trajectory.json"` file for each individual run.
///
void parse(
    InputParser<monte::jsonResultsIO<config_type>> &parser,
    monte::StateSamplingFunctionMap<config_type> const &sampling_functions,
    monte::ResultsAnalysisFunctionMap<config_type> const &analysis_functions) {
  std::string output_dir = "output";
  parser.optional(output_dir, "output_dir");

  bool write_observations = false;
  parser.optional(write_observations, "write_observations");

  bool write_trajectory = false;
  parser.optional(write_trajectory, "write_trajectory");

  if (parser.valid()) {
    parser.value = std::make_unique<monte::jsonResultsIO<config_type>>(
        output_dir, sampling_functions, analysis_functions, write_trajectory,
        write_observations);
  }
}

}  // namespace clexmonte
}  // namespace CASM
