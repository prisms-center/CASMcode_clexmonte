#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/functions.hh"
#include "casm/clexmonte/canonical/run.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "gtest/gtest.h"
#include "misc.hh"
#include "testdir.hh"

using namespace CASM;

TEST(canonical_run_test, Test1) {
  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path output_dir_relpath = "output";

  fs::path test_dir =
      fs::current_path() / "CASM_test_projects" / "canonical_run_test";
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

  /// Make sampling & analysis functions
  auto sampling_functions =
      clexmonte::canonical::make_sampling_functions(system);
  auto analysis_functions =
      clexmonte::canonical::make_analysis_functions(system);

  /// Parse and construct RunParams
  auto run_subparser = parser.subparse<clexmonte::RunParams>(
      "kwargs", system, sampling_functions, analysis_functions);
  std::runtime_error run_error_if_invalid{
      "Error reading canonical Monte Carlo run JSON input"};
  report_and_throw_if_invalid(*run_subparser, CASM::log(),
                              run_error_if_invalid);
  EXPECT_TRUE(run_subparser->valid());

  clexmonte::RunParams &run_params = *run_subparser->value;

  // ### Random number generator engine pointer (empty - will be seeded by
  // random device)
  std::shared_ptr<std::mt19937_64> random_number_engine;

  clexmonte::canonical::run(
      system, run_params.sampling_functions, run_params.analysis_functions,
      *run_params.state_generator, run_params.sampling_params,
      run_params.completion_check_params, *run_params.results_io,
      random_number_engine);

  EXPECT_TRUE(fs::exists(test_dir / "output"));
  EXPECT_TRUE(fs::exists(test_dir / "output" / "summary.json"));
  EXPECT_EQ(test::file_count(test_dir / "output"), 1);

  fs::remove(test_dir / "output" / "summary.json");
  fs::remove(test_dir / "output");
}
