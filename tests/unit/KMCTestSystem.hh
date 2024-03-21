#ifndef CASM_unittest_KMCTestSystem
#define CASM_unittest_KMCTestSystem

#include <filesystem>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/global/definitions.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

namespace test {

using namespace CASM;

/// NOTE:
/// - This test fixture is designed to copy data to the same directory
///   each time, so that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
/// CASM_test_projects/<_test_dir_name>
class KMCTestSystem : public testing::Test {
 protected:
  std::string project_name;
  fs::path test_data_dir;
  fs::path test_dir;
  fs::copy_options copy_options;
  jsonParser json;

  std::shared_ptr<clexmonte::System> system;

  /// \brief Default test project - if you use this, don't copy other
  ///     files into the test project, just use the default test fixture
  KMCTestSystem();

  /// \brief Constructor
  KMCTestSystem(std::string _project_name, std::string _test_dir_name,
                fs::path _input_file_path);

  void setup_input_files(bool use_sparse_format_eci);

  /// \brief Copy formation_energy Clexulator and ECI to test_dir and update
  ///     input json with location
  void set_clex(std::string clex_name, std::string bset_name,
                fs::path eci_relpath);

  void copy_local_clexulator(fs::path src_basis_sets_dir,
                             fs::path dest_basis_sets_dir,
                             std::string bset_name,
                             std::string clexulator_basename);

  /// \brief Copy local clexulator and eci to test_dir and update input json
  /// with location
  void set_local_basis_set(std::string bset_name);

  void set_event(std::string event_name, std::string kra_eci_relpath,
                 std::string freq_eci_relpath);

  void write_input();

  void make_system();
};

}  // namespace test

#endif
