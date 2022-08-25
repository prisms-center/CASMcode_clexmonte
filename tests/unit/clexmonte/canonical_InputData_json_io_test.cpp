#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/io/InputData.hh"
#include "casm/clexmonte/canonical/io/json/InputData_json_io.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace CASM;

TEST(canonical_InputData_json_io_test, Test1) {
  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path output_dir_relpath = "output";

  test::TmpDir tmp_dir;
  fs::create_directories(tmp_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                tmp_dir / clexulator_src_relpath);
  fs::copy_file(test_data_dir / eci_relpath, tmp_dir / eci_relpath);
  fs::create_directories(tmp_dir / output_dir_relpath);

  jsonParser json(test_data_dir / "input_1.json");
  json["kwargs"]["system"]["basis_sets"]["formation_energy"]["source"] =
      (tmp_dir / clexulator_src_relpath).string();
  json["kwargs"]["system"]["clex"]["formation_energy"]["coefficients"] =
      (tmp_dir / eci_relpath).string();
  json["kwargs"]["results_io"]["kwargs"]["output_dir"] =
      (tmp_dir / output_dir_relpath).string();

  ParentInputParser parser(json);
  auto subparser = parser.subparse<clexmonte::canonical::InputData>("kwargs");
  std::runtime_error error_if_invalid{
      "Error reading canonical Monte Carlo JSON input"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  EXPECT_TRUE(subparser->valid());
}
