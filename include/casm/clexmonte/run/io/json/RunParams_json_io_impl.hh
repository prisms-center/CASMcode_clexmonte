#ifndef CASM_clexmonte_run_RunParams_json_io_impl
#define CASM_clexmonte_run_RunParams_json_io_impl

#include "casm/clexmonte/misc/subparse_from_file.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/run/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"
#include "casm/monte/state/StateGenerator.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse canonical Monte Carlo input file
///
/// Input file summary:
/// \code
/// {
///     "state_generation": <monte::StateGenerator>
///         Specifies a "path" of input states at which to run Monte Carlo
///         calculations. Each state is an initial configuration and set of
///         thermodynamic conditions (temperature, chemical potential,
///         composition, etc.).
///     "random_number_generator":
///         (Future) Options controlling the random number generator.
///     "sampling_fixtures": JSON object
///         A JSON object, whose keys are labels and values are paths to
///         input files for sampling fixtures. A Monte Carlo run continues
///         until all sampling fixtures are completed.
///     "before_first_run": optional JSON object = null
///         An optional JSON object, whose keys are labels and values are
///         paths to input files for sampling fixtures. If included, the
///         requested run will be performed at the initial conditions as
///         a preliminary step before the actual first run begins. This
///         may be useful when not running in automatic convergence mode.
///     "before_each_run": optional JSON object = null
///         An optional JSON object, whose keys are labels and values are
///         paths to input files for sampling fixtures. If included, the
///         requested run will be performed as a preliminary step before
///         each actual run begins. This may be useful when not running
///         in automatic convergence mode.
///     "completed_runs/save_all_initial_states": bool = false
///         If true, save initial states for analysis, output, or state
///         generation.
///     "completed_runs/save_all_final_states": bool = false
///         If true, save final states for analysis, output, or state
///         generation.
///     "completed_runs/save_last_final_state": bool = true
///         If true, save final state for last run for analysis, output
///         or state generation.
///     "completed_runs/write_initial_states": bool = false
///         If true, write saved initial states to completed_runs.json.
///     "completed_runs/write_final_states": bool = false
///         If true, write saved final states to completed_runs.json.
///     "completed_runs/output_dir": str = ""
///         If not empty, name of a directory in which to write
///         completed_runs.json. If empty, completed_runs.json is not
///         written, which means restarts are not possible.
///     "global_cutoff": bool = true
///         If true, the entire run is stopped when any sampling fixture
///         is completed. Otherwise, all fixtures must complete for the
///         run to be completed.
/// }
/// \endcode
template <typename EngineType>
void parse(InputParser<RunParams<EngineType>> &parser,
           std::vector<fs::path> search_path,
           std::shared_ptr<EngineType> engine,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions,
           MethodParserMap<state_generator_type> const &state_generator_methods,
           MethodParserMap<results_io_type> const &results_io_methods,
           bool time_sampling_allowed) {
  /// TODO: "random_number_generator":
  ///     (Future) Options controlling the random number generator.

  // Construct state generator
  auto state_generator_subparser =
      parser.template subparse<state_generator_type>("state_generation",
                                                     state_generator_methods);

  // Construct sampling fixture parameters
  auto _parse_sampling_fixtures = [&](std::string key, bool is_required) {
    std::vector<sampling_fixture_params_type> sampling_fixture_params;
    if (parser.self.contains(key) && !parser.self[key].is_null()) {
      auto it = parser.self[key].begin();
      auto end = parser.self[key].end();
      for (; it != end; ++it) {
        std::string label = it.name();
        std::shared_ptr<InputParser<sampling_fixture_params_type>> subparser;
        if (it->is_obj()) {
          subparser = parser.template subparse<sampling_fixture_params_type>(
              fs::path(key) / label, label, sampling_functions,
              analysis_functions, results_io_methods, time_sampling_allowed);
        } else if (it->is_string()) {
          subparser = subparse_from_file<sampling_fixture_params_type>(
              parser, fs::path(key) / label, search_path, label,
              sampling_functions, analysis_functions, results_io_methods,
              time_sampling_allowed);
        } else {
          parser.insert_error(fs::path(key) / label,
                              "Error: must be a file name or JSON object");
          continue;
        }
        if (subparser->valid()) {
          sampling_fixture_params.push_back(*subparser->value);
        }
      }
    } else if (is_required) {
      std::stringstream msg;
      msg << "Error: '" << key << "' is required.";
      parser.insert_error(key, msg.str());
    }
    return sampling_fixture_params;
  };

  bool is_required;
  std::vector<sampling_fixture_params_type> sampling_fixture_params =
      _parse_sampling_fixtures("sampling_fixtures", is_required = true);
  std::vector<sampling_fixture_params_type> before_first_run =
      _parse_sampling_fixtures("before_first_run", is_required = false);
  std::vector<sampling_fixture_params_type> before_each_run =
      _parse_sampling_fixtures("before_each_run", is_required = false);

  monte::RunManagerParams run_manager_params;
  parser.optional_else(run_manager_params.do_save_all_initial_states,
                       "completed_runs/save_all_initial_states", false);
  parser.optional_else(run_manager_params.do_save_all_final_states,
                       "completed_runs/save_all_final_states", false);
  parser.optional_else(run_manager_params.do_save_all_final_states,
                       "completed_runs/save_all_final_states", false);
  parser.optional_else(run_manager_params.do_save_last_final_state,
                       "completed_runs/save_last_final_state", true);
  parser.optional_else(run_manager_params.do_write_initial_states,
                       "completed_runs/write_initial_states", false);
  parser.optional_else(run_manager_params.do_write_final_states,
                       "completed_runs/write_final_states", false);
  std::string output_dir;
  parser.optional(output_dir, "completed_runs/output_dir");
  run_manager_params.output_dir = fs::path(output_dir);
  parser.optional_else(run_manager_params.global_cutoff, "global_cutoff", true);

  if (parser.valid()) {
    parser.value = std::make_unique<RunParams<EngineType>>(
        engine, std::move(state_generator_subparser->value), run_manager_params,
        sampling_fixture_params, before_first_run, before_each_run);
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
