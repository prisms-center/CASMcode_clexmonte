#include "ZrOTestSystem.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/global/filesystem.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

namespace test {

/// \brief Default test project - if you use this, don't copy other
///     files into the test project, just use the default test fixture
ZrOTestSystem::ZrOTestSystem()
    : ZrOTestSystem(
          "ZrOTestSystem_default",
          test::data_dir("clexmonte") / "ZrOTestSystem" / "system.json") {
  std::string clex_name = "formation_energy";
  std::string bset_name = "formation_energy";
  fs::path eci_relpath = "formation_energy_eci.json";
  set_clex(clex_name, bset_name, eci_relpath);
  make_system();
}

/// \brief Use this constructor if you want to test additional clex, etc.
ZrOTestSystem::ZrOTestSystem(std::string _test_dir_name,
                             fs::path _input_file_path)
    : project_name("ZrO"),
      test_data_dir(test::data_dir("clexmonte") / "ZrOTestSystem"),
      test_dir(fs::current_path() / "CASM_test_projects" / _test_dir_name),
      copy_options(fs::copy_options::skip_existing),
      system_json(_input_file_path) {
  fs::create_directories(test_dir);
}

/// \brief Copy formation_energy Clexulator and ECI to test_dir and update
///     input json with location
///
/// Notes:
/// - Assumes Clexulator file is:
///   - basis_sets/bset.<bset_name>/<project_name>_Clexulator_<bset_name>.cc
/// - ECI files is <eci_relpath>
void ZrOTestSystem::set_clex(std::string clex_name, std::string bset_name,
                             fs::path eci_relpath) {
  fs::path clexulator_src_relpath =
      fs::path("basis_sets") / ("bset." + bset_name) /
      (project_name + "_Clexulator_" + bset_name + ".cc");

  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::create_directories(test_dir / eci_relpath.parent_path());
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);

  system_json["basis_sets"][bset_name]["source"] =
      (test_dir / clexulator_src_relpath).string();
  system_json["clex"][clex_name]["basis_set"] = bset_name;
  system_json["clex"][clex_name]["coefficients"] =
      (test_dir / eci_relpath).string();
}

void ZrOTestSystem::make_system() {
  std::vector<fs::path> search_path;
  InputParser<clexmonte::System> parser(system_json, search_path);
  std::runtime_error error_if_invalid{"Error reading ZrOTestSystem data"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  system = std::shared_ptr<clexmonte::System>(std::move(parser.value));
}

}  // namespace test
