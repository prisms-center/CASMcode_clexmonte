#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/FixedConfigGenerator.hh"
#include "casm/clexmonte/run/StateGenerator.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/ConfigGenerator_json_io.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io_impl.hh"
#include "casm/clexmonte/run/io/json/StateGenerator_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/run_management/io/ResultsIO.hh"
#include "casm/monte/run_management/io/json/ResultsIO_json_io.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/monte/run_management/io/json/jsonResultsIO_impl.hh"

namespace CASM {
namespace clexmonte {

MethodParserMap<config_generator_type> standard_config_generator_methods(
    std::shared_ptr<system_type> const &system) {
  MethodParserFactory<config_generator_type> cf;
  MethodParserMap<config_generator_type> config_generator_methods;
  config_generator_methods.insert(
      cf.make<FixedConfigGenerator>("fixed", system)
      // To add additional config generators:
      // cf.make<DerivedClassName>("<name>", ...args...),
  );
  return config_generator_methods;
}

MethodParserMap<results_io_type> standard_results_io_methods() {
  MethodParserFactory<results_io_type> f;
  MethodParserMap<results_io_type> results_io_methods;

  results_io_methods.insert(
      f.template make<monte::jsonResultsIO<results_type>>("json")
      // To add additional state generators:
      // f.make<DerivedClassName>("<name>", ...args...),
  );
  return results_io_methods;
}

}  // namespace clexmonte
}  // namespace CASM
