#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/run/io/RunParams.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io_impl.hh"
#include "casm/clexmonte/semigrand_canonical/calculator.hh"
#include "casm/clexmonte/semigrand_canonical/json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "gtest/gtest.h"
#include "misc.hh"
#include "testdir.hh"

using namespace CASM;

/// Simple test of semi-grand canonical Monte Carlo
TEST(semigrand_canonical_run_test, Test1) {
  std::vector<fs::path> search_path;

  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path output_dir_relpath = "output";

  fs::path test_dir = fs::current_path() / "CASM_test_projects" /
                      "semigrand_canonical_run_test";
  fs::copy_options copy_options = fs::copy_options::skip_existing;
  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);
  fs::copy_file(test_data_dir / "system.json", test_dir / eci_relpath,
                copy_options);

  /// Parse and construct system
  jsonParser system_json(test_data_dir / "system.json");
  system_json["basis_sets"]["formation_energy"]["source"] =
      (test_dir / clexulator_src_relpath).string();
  system_json["clex"]["formation_energy"]["coefficients"] =
      (test_dir / eci_relpath).string();
  InputParser<clexmonte::System> system_parser(system_json, search_path);
  std::runtime_error system_error_if_invalid{
      "Error reading semi-grand canonical Monte Carlo system JSON input"};
  report_and_throw_if_invalid(system_parser, CASM::log(),
                              system_error_if_invalid);

  std::shared_ptr<clexmonte::System> system(system_parser.value.release());

  // Make calculation object:
  typedef clexmonte::semigrand_canonical::SemiGrandCanonical_mt19937_64
      calculation_type;
  typedef calculation_type::engine_type engine_type;
  auto calculation = std::make_shared<calculation_type>(system);

  // Make random number engine
  std::shared_ptr<engine_type> engine = std::make_shared<engine_type>();

  /// Make state sampling & analysis functions
  auto sampling_functions =
      calculation_type::standard_sampling_functions(calculation);
  auto json_sampling_functions =
      calculation_type::standard_json_sampling_functions(calculation);
  auto analysis_functions =
      calculation_type::standard_analysis_functions(calculation);
  auto modifying_functions =
      calculation_type::standard_modifying_functions(calculation);

  /// Make config generator / state generator / results_io JSON parsers
  clexmonte::semigrand_canonical::SemiGrandCanonicalConditions const
      *conditions_ptr;
  auto config_generator_methods =
      clexmonte::standard_config_generator_methods(calculation->system);
  auto state_generator_methods = clexmonte::standard_state_generator_methods(
      calculation->system, modifying_functions, config_generator_methods,
      conditions_ptr);
  auto results_io_methods = clexmonte::standard_results_io_methods();

  /// Parse and construct run parameters
  std::cout << "here 0" << std::endl;
  jsonParser run_params_json(test_data_dir / "run_params_sgc_complete.json");
  run_params_json["sampling_fixtures"]["thermo"]["results_io"]["kwargs"]
                 ["output_dir"] =
                     (test_dir / output_dir_relpath / "thermo").string();
  std::cout << "here 1" << std::endl;
  InputParser<clexmonte::RunParams<std::mt19937_64>> run_params_parser(
      run_params_json, search_path, engine, sampling_functions,
      json_sampling_functions, analysis_functions, state_generator_methods,
      results_io_methods, calculation->time_sampling_allowed, conditions_ptr);
  std::cout << "here 2" << std::endl;
  std::runtime_error run_params_error_if_invalid{
      "Error reading Monte Carlo run parameters JSON input"};
  report_and_throw_if_invalid(run_params_parser, CASM::log(),
                              run_params_error_if_invalid);

  std::cout << "here 3" << std::endl;
  clexmonte::RunParams<std::mt19937_64> &run_params = *run_params_parser.value;

  std::cout << "begin run_series" << std::endl;
  clexmonte::run_series(*calculation, engine, *run_params.state_generator,
                        run_params.sampling_fixture_params,
                        run_params.global_cutoff, run_params.before_first_run,
                        run_params.before_each_run);

  EXPECT_TRUE(fs::exists(test_dir / "output"));
  EXPECT_TRUE(fs::exists(test_dir / "output" / "thermo"));
  EXPECT_TRUE(fs::exists(test_dir / "output" / "thermo" / "summary.json"));
  EXPECT_EQ(test::file_count(test_dir / "output" / "thermo"), 1);

  fs::remove(test_dir / "output" / "thermo" / "summary.json");
  fs::remove(test_dir / "output" / "thermo");
  fs::remove(test_dir / "output");
}

/// Test semi-grand canonical Monte Carlo with 2 sampling fixtures, specified
/// using seperate files.
/// - First sampling fixture, labeled "thermo_period1":  sampling every
///   1 pass and specified using file "thermo_sampling.period1.json"
/// - Second sampling fixture, labeled "thermo_period10": sampling every
///   10 passes and specified using file "thermo_sampling.period10.json"
TEST(semigrand_canonical_run_test, Test2) {
  std::vector<fs::path> search_path;

  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path output_dir_relpath = "output";

  fs::path test_dir = fs::current_path() / "CASM_test_projects" /
                      "semigrand_canonical_run_test";
  fs::copy_options copy_options = fs::copy_options::skip_existing;
  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);
  fs::copy_file(test_data_dir / "system.json", test_dir / eci_relpath,
                copy_options);

  /// Construct "thermo_sampling.period1.json
  jsonParser p1_json(test_data_dir / "thermo_sampling.period1.json");
  p1_json["results_io"]["kwargs"]["output_dir"] =
      (test_dir / output_dir_relpath / "thermo_period1").string();
  p1_json.write(test_dir / "thermo_sampling.period1.json");

  // for debugging:
  // p1_json["log"]["file"] =
  //     (test_dir / output_dir_relpath / "thermo_period1" / "status.json")
  //         .string();

  /// Construct "thermo_sampling.period10.json
  jsonParser p10_json(test_data_dir / "thermo_sampling.period10.json");
  p10_json["results_io"]["kwargs"]["output_dir"] =
      (test_dir / output_dir_relpath / "thermo_period10").string();
  p10_json.write(test_dir / "thermo_sampling.period10.json");

  // for debugging:
  // p10_json["log"]["file"] =
  //     (test_dir / output_dir_relpath / "thermo_period10" / "status.json")
  //         .string();

  /// Parse and construct system
  jsonParser system_json(test_data_dir / "system.json");
  system_json["basis_sets"]["formation_energy"]["source"] =
      (test_dir / clexulator_src_relpath).string();
  system_json["clex"]["formation_energy"]["coefficients"] =
      (test_dir / eci_relpath).string();
  InputParser<clexmonte::System> system_parser(system_json, search_path);
  std::runtime_error system_error_if_invalid{
      "Error reading semi-grand canonical Monte Carlo system JSON input"};
  report_and_throw_if_invalid(system_parser, CASM::log(),
                              system_error_if_invalid);

  std::shared_ptr<clexmonte::System> system(system_parser.value.release());

  // Make calculation object:
  typedef clexmonte::semigrand_canonical::SemiGrandCanonical_mt19937_64
      calculation_type;
  typedef calculation_type::engine_type engine_type;
  auto calculation = std::make_shared<calculation_type>(system);

  // Make random number engine
  std::shared_ptr<engine_type> engine = std::make_shared<engine_type>();

  /// Make state sampling & analysis functions
  auto sampling_functions =
      calculation_type::standard_sampling_functions(calculation);
  auto json_sampling_functions =
      calculation_type::standard_json_sampling_functions(calculation);
  auto analysis_functions =
      calculation_type::standard_analysis_functions(calculation);
  auto modifying_functions =
      calculation_type::standard_modifying_functions(calculation);

  /// Make config generator / state generator / results_io JSON parsers
  clexmonte::semigrand_canonical::SemiGrandCanonicalConditions const
      *conditions_ptr = nullptr;
  auto config_generator_methods =
      clexmonte::standard_config_generator_methods(calculation->system);
  auto state_generator_methods = clexmonte::standard_state_generator_methods(
      calculation->system, modifying_functions, config_generator_methods,
      conditions_ptr);
  auto results_io_methods = clexmonte::standard_results_io_methods();

  /// Parse and construct run parameters
  jsonParser run_params_json(test_data_dir / "run_params_sgc_by_file.json");
  run_params_json["sampling_fixtures"]["thermo_period1"] =
      (test_dir / "thermo_sampling.period1.json").string();
  run_params_json["sampling_fixtures"]["thermo_period10"] =
      (test_dir / "thermo_sampling.period10.json").string();
  InputParser<clexmonte::RunParams<std::mt19937_64>> run_params_parser(
      run_params_json, search_path, engine, sampling_functions,
      json_sampling_functions, analysis_functions, state_generator_methods,
      results_io_methods, calculation->time_sampling_allowed, conditions_ptr);
  std::runtime_error run_params_error_if_invalid{
      "Error reading Monte Carlo run parameters JSON input"};
  report_and_throw_if_invalid(run_params_parser, CASM::log(),
                              run_params_error_if_invalid);

  clexmonte::RunParams<std::mt19937_64> &run_params = *run_params_parser.value;

  clexmonte::run_series(*calculation, engine, *run_params.state_generator,
                        run_params.sampling_fixture_params,
                        run_params.global_cutoff, run_params.before_first_run,
                        run_params.before_each_run);

  EXPECT_TRUE(fs::exists(test_dir / "output"));

  EXPECT_TRUE(fs::exists(test_dir / "output" / "thermo_period1"));
  EXPECT_TRUE(
      fs::exists(test_dir / "output" / "thermo_period1" / "summary.json"));
  EXPECT_EQ(test::file_count(test_dir / "output" / "thermo_period1"), 1);
  fs::remove(test_dir / "output" / "thermo_period1" / "summary.json");
  fs::remove(test_dir / "output" / "thermo_period1");

  EXPECT_TRUE(fs::exists(test_dir / "output" / "thermo_period10"));
  EXPECT_TRUE(
      fs::exists(test_dir / "output" / "thermo_period10" / "summary.json"));
  EXPECT_EQ(test::file_count(test_dir / "output" / "thermo_period10"), 1);
  fs::remove(test_dir / "output" / "thermo_period10" / "summary.json");
  fs::remove(test_dir / "output" / "thermo_period10");

  fs::remove(test_dir / "output");

  fs::remove(test_dir / "thermo_sampling.period1.json");
  fs::remove(test_dir / "thermo_sampling.period10.json");
}
