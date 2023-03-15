#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/run/io/json/ResultsIO_json_io.hh"
#include "casm/clexmonte/run/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/monte/results/io/ResultsIO.hh"
#include "casm/monte/results/io/json/jsonResultsIO_impl.hh"
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
      cf.make<monte::FixedConfigGenerator<config_type>>("fixed", system)
      // To add additional config generators:
      // cf.make<DerivedClassName>("<name>", ...args...),
  );
  return config_generator_methods;
}

MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    std::map<std::string, state_modifying_function_type> const
        &modifying_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods) {
  MethodParserFactory<state_generator_type> sf;
  MethodParserMap<state_generator_type> state_generator_methods;
  state_generator_methods.insert(
      sf.make<monte::IncrementalConditionsStateGenerator<Configuration>>(
          "incremental", system, modifying_functions, config_generator_methods)
      // To add additional state generators:
      // sf.make<DerivedClassName>("<name>", ...args...),
  );
  return state_generator_methods;
}

MethodParserMap<results_io_type> standard_results_io_methods(
    std::map<std::string, state_sampling_function_type> const
        &sampling_functions,
    std::map<std::string, results_analysis_function_type> const
        &analysis_functions) {
  MethodParserFactory<results_io_type> f;
  MethodParserMap<results_io_type> results_io_methods;

  results_io_methods.insert(f.template make<monte::jsonResultsIO<results_type>>(
      "json", sampling_functions, analysis_functions)
                            // To add additional state generators:
                            // f.make<DerivedClassName>("<name>", ...args...),
  );
  return results_io_methods;
}

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
///     "random_number_generator": <monte::RandomNumberGenerator>
///         (Future) Options controlling the random number generator.
///     "sampling_fixtures": JSON object
///         A JSON object, whose keys are labels and values are paths to
///         input files for sampling fixtures. A Monte Carlo run continues
///         until all sampling fixtures are completed.
///     "save_all_initial_states": bool = false
///         If true, save initial states for analysis, output, or state
///         generation.
///     "save_all_final_states": bool = false
///         If true, save final states for analysis, output, or state
///         generation.
///     "save_last_final_state": bool = true
///         If true, save final state for last run for analysis, output
///         or state generation.
///     "write_initial_states": bool = false
///         If true, write saved initial states to completed_runs.json.
///     "write_final_states": bool = false
///         If true, write saved final states to completed_runs.json.
///     "output_dir": str = ""
///         If not empty, name of a directory in which to write
///         completed_runs.json. If empty, completed_runs.json is not
///         written, which means restarts are not possible.
///
/// }
/// \endcode
void parse(InputParser<RunParams> &parser,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions,
           MethodParserMap<state_generator_type> const &state_generator_methods,
           MethodParserMap<results_io_type> const &results_io_methods) {
  // Construct state generator
  auto state_generator_subparser = parser.subparse<state_generator_type>(
      "state_generation", state_generator_methods);

  // Construct sampling fixture parameters
  std::vector<sampling_fixture_params_type> sampling_fixture_params;
  if (parser.self.contains("sampling_fixtures")) {
    auto it = parser.self["sampling_fixtures"].begin();
    auto end = parser.self["sampling_fixtures"].end();
    for (; it != end; ++it) {
      std::string label = it.name();
      std::shared_ptr<InputParser<sampling_fixture_params_type>> subparser;
      if (it->is_obj()) {
        subparser = parser.subparse<sampling_fixture_params_type>(
            fs::path("sampling_fixtures") / label, label, sampling_functions,
            analysis_functions, results_io_methods);
      } else if (it->is_string()) {
        subparser = parser.subparse_from_file<sampling_fixture_params_type>(
            fs::path("sampling_fixtures") / label, label, sampling_functions,
            analysis_functions, results_io_methods);
      } else {
        parser.insert_error(fs::path("sampling_params") / label,
                            "Error: must be a file name or JSON object");
        continue;
      }
      if (subparser->valid()) {
        sampling_fixture_params.push_back(*subparser->value);
      }
    }
  } else {
    parser.insert_error("sampling_fixtures", "Error: no 'sampling_fixtures'");
  }

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

  if (parser.valid()) {
    parser.value = std::make_unique<RunParams>(
        std::move(state_generator_subparser->value), run_manager_params,
        sampling_fixture_params);
  }
}

}  // namespace clexmonte
}  // namespace CASM
