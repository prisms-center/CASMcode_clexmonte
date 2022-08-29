#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/io/InputData.hh"
#include "casm/clexmonte/canonical/io/json/InputData_json_io.hh"
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
  auto subparser = parser.subparse<clexmonte::canonical::InputData>("kwargs");
  std::runtime_error error_if_invalid{
      "Error reading canonical Monte Carlo JSON input"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_TRUE(subparser->valid());

  clexmonte::canonical::InputData &input_data = *subparser->value;
  clexmonte::canonical::run(input_data);

  EXPECT_TRUE(fs::exists(test_dir / "output"));
  EXPECT_TRUE(fs::exists(test_dir / "output" / "summary.json"));
  EXPECT_EQ(test::file_count(test_dir / "output"), 1);

  fs::remove(test_dir / "output" / "summary.json");
  fs::remove(test_dir / "output");
}
