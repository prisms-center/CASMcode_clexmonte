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

}  // namespace clexmonte
}  // namespace CASM
