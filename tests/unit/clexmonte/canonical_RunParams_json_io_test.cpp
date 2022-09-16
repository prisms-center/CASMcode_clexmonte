#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace CASM;

TEST(canonical_RunParams_json_io_test, Test1) {
  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path output_dir_relpath = "output";

  fs::path test_dir = fs::current_path() / "CASM_test_projects" /
                      "canonical_RunParams_json_io_test";
  fs::copy_options copy_options = fs::copy_options::skip_existing;
  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);
  fs::create_directories(test_dir / output_dir_relpath);

  jsonParser json(test_data_dir / "input_1.json");
  json["kwargs"]["system"]["basis_sets"]["formation_energy"]["source"] =
      (test_dir / clexulator_src_relpath).string();
  json["kwargs"]["system"]["clex"]["formation_energy"]["coefficients"] =
      (test_dir / eci_relpath).string();
  json["kwargs"]["results_io"]["kwargs"]["output_dir"] =
      (test_dir / output_dir_relpath).string();

  ParentInputParser parser(json);

  /// Parse and construct system
  auto system_subparser =
      parser.subparse<clexmonte::System>(fs::path("kwargs") / "system");
  std::runtime_error system_error_if_invalid{
      "Error reading canonical Monte Carlo system JSON input"};
  report_and_throw_if_invalid(*system_subparser, CASM::log(),
                              system_error_if_invalid);

  std::shared_ptr<clexmonte::System> system(system_subparser->value.release());

  // Make calculation object:
  typedef clexmonte::canonical::Canonical_mt19937_64 calculation_type;
  auto calculation = std::make_shared<calculation_type>(system);

  /// Make state sampling & analysis functions
  auto sampling_functions =
      calculation_type::standard_sampling_functions(calculation);
  auto analysis_functions =
      calculation_type::standard_analysis_functions(calculation);

  /// Make config generator / state generator / results_io JSON parsers
  auto config_generator_methods =
      clexmonte::standard_config_generator_methods(calculation->system);
  auto state_generator_methods = clexmonte::standard_state_generator_methods(
      calculation->system, sampling_functions, config_generator_methods);
  auto results_io_methods = clexmonte::standard_results_io_methods(
      sampling_functions, analysis_functions);

  /// Parse and construct run parameters
  auto run_params_subparser = parser.subparse<clexmonte::RunParams>(
      fs::path("kwargs"), system, sampling_functions, analysis_functions,
      state_generator_methods, results_io_methods);
  std::runtime_error run_params_error_if_invalid{
      "Error reading Monte Carlo run parameters JSON input"};
  report_and_throw_if_invalid(*run_params_subparser, CASM::log(),
                              run_params_error_if_invalid);

  EXPECT_TRUE(run_params_subparser->valid());

  // -- remove "output" dir
  fs::remove(test_dir / "output");
}
