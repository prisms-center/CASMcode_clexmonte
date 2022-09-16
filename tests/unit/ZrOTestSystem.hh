#ifndef CASM_unittest_ZrOTestSystem
#define CASM_unittest_ZrOTestSystem

#include <filesystem>

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/global/filesystem.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

namespace test {
using namespace CASM;
using namespace CASM::monte;
using namespace CASM::clexmonte;

/// configure with (in order of priority):
/// - CASM_CXX, CXX, default="c++"
/// - CASM_CXXFLAGS, default="-O3 -Wall -fPIC --std=c++17"
/// - CASM_SOFLAGS, default="-shared"
class ZrOTestSystem : public testing::Test {
 public:
  std::string project_name;
  fs::path test_data_dir;
  fs::path test_dir;
  fs::copy_options copy_options;
  jsonParser system_json;

  std::shared_ptr<clexmonte::System> system;

  /// \brief Default test project - if you use this, don't copy other
  ///     files into the test project, just use the default test fixture
  ZrOTestSystem();

  /// \brief Use this constructor if you want to test additional clex, etc.
  ZrOTestSystem(std::string _test_dir_name, fs::path _input_file_path);

  /// \brief Copy formation_energy Clexulator and ECI to test_dir and update
  ///     input json with location
  void set_clex(std::string clex_name, std::string bset_name,
                fs::path eci_relpath);

  void make_system();
};

}  // namespace test

#endif
