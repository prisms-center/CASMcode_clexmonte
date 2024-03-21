#ifndef CASM_clexmonte_run_RunParams_json_io_impl
#define CASM_clexmonte_run_RunParams_json_io_impl

#include "casm/clexmonte/misc/subparse_from_file.hh"
#include "casm/clexmonte/run/FixedConfigGenerator.hh"
#include "casm/clexmonte/run/IncrementalConditionsStateGenerator.hh"
#include "casm/clexmonte/run/StateGenerator.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io_impl.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"

namespace CASM {
namespace clexmonte {

template <typename ConditionsType>
MethodParserMap<state_generator_type> standard_state_generator_methods(
    std::shared_ptr<system_type> const &system,
    StateModifyingFunctionMap const &modifying_functions,
    MethodParserMap<config_generator_type> const &config_generator_methods,
    ConditionsType const *ptr) {
  MethodParserFactory<state_generator_type> sf;
  MethodParserMap<state_generator_type> state_generator_methods;
  state_generator_methods.insert(
      sf.make<IncrementalConditionsStateGenerator>(
          "incremental", system, modifying_functions, config_generator_methods,
          ptr)
      // To add additional state generators:
      // sf.make<DerivedClassName>("<name>", ...args...),
  );
  return state_generator_methods;
}

/// \brief Parse canonical Monte Carlo input file
///
/// Input file summary:
/// \code
/// {
///     "state_generation": <StateGenerator>
///         Specifies a "path" of input states at which to run Monte Carlo
///         calculations. Each state is an initial configuration and set of
///         thermodynamic conditions (temperature, chemical potential,
///         composition, etc.).
///     "random_number_generator":
///         (TODO) Options controlling the random number generator.
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
///     "global_cutoff": bool = true
///         If true, the entire run is stopped when any sampling fixture
///         is completed. Otherwise, all fixtures must complete for the
///         run to be completed.
/// }
/// \endcode
template <typename EngineType, typename ConditionsType>
void parse(InputParser<RunParams<EngineType>> &parser,
           std::vector<fs::path> search_path,
           std::shared_ptr<EngineType> engine,
           monte::StateSamplingFunctionMap const &sampling_functions,
           monte::jsonStateSamplingFunctionMap const &json_sampling_functions,
           monte::ResultsAnalysisFunctionMap<config_type, statistics_type> const
               &analysis_functions,
           MethodParserMap<state_generator_type> const &state_generator_methods,
           MethodParserMap<results_io_type> const &results_io_methods,
           bool time_sampling_allowed, ConditionsType const *ptr) {
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
              json_sampling_functions, analysis_functions, results_io_methods,
              time_sampling_allowed);
        } else if (it->is_string()) {
          subparser = subparse_from_file<sampling_fixture_params_type>(
              parser, fs::path(key) / label, search_path, label,
              sampling_functions, json_sampling_functions, analysis_functions,
              results_io_methods, time_sampling_allowed);
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

  bool global_cutoff;
  parser.optional_else(global_cutoff, "global_cutoff", true);

  if (parser.valid()) {
    parser.value = std::make_unique<RunParams<EngineType>>(
        engine, std::move(state_generator_subparser->value),
        sampling_fixture_params, global_cutoff, before_first_run,
        before_each_run);
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
