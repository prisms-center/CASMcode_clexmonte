#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/definitions.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

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

  /// \param _project_name Name of project test files to use
  /// (clexmonte/data/<_project_name>) \param _test_dir_name Name of directory
  /// where test files should be copied and tested
  /// (CASM_test_projects/<_test_dir_name>) \param _input_file_path Path to
  /// template input file that will be updated with the
  ///     paths to the copied Clexulator and ECI files
  ///     (test::data_dir("clexmonte") / "kmc" / "system_1.json")
  KMCTestSystem(std::string _project_name, std::string _test_dir_name,
                fs::path _input_file_path)
      : project_name(_project_name),
        test_data_dir(test::data_dir("clexmonte") / _project_name),
        test_dir(fs::current_path() / "CASM_test_projects" / _test_dir_name),
        copy_options(fs::copy_options::skip_existing),
        json(_input_file_path) {
    std::cout << "KMCTestSystem::test_dir: " << test_dir << std::endl;
  }

  /// \brief Copy formation_energy Clexulator and ECI to test_dir and update
  ///     input json with location
  ///
  /// Notes:
  /// - Assumes Clexulator file is:
  ///   - basis_sets/bset.<bset_name>/<project_name>_Clexulator_<bset_name>.cc
  /// - ECI files is <eci_relpath>
  void set_clex(std::string clex_name, std::string bset_name,
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

    json["kwargs"]["system"][clex_name]["source"] =
        (test_dir / clexulator_src_relpath).string();
    json["kwargs"]["system"][clex_name]["coefficients"] =
        (test_dir / eci_relpath).string();
  }

  void copy_local_clexulator(fs::path src_basis_sets_dir,
                             fs::path dest_basis_sets_dir,
                             std::string bset_name,
                             std::string clexulator_basename) {
    fs::path src_dir = src_basis_sets_dir / (std::string("bset.") + bset_name);
    fs::path dest_dir =
        dest_basis_sets_dir / (std::string("bset.") + bset_name);

    // equivalents
    Index i = 0;
    fs::path equiv_dir = fs::path(std::to_string(i));
    while (fs::exists(src_dir / equiv_dir)) {
      std::string src_filename = clexulator_basename + "_" + bset_name + "_" +
                                 std::to_string(i) + ".cc";
      if (!fs::exists(src_dir / equiv_dir / src_filename)) {
        break;
      }
      fs::create_directories(dest_dir / equiv_dir);
      fs::copy_file(src_dir / equiv_dir / src_filename,
                    dest_dir / equiv_dir / src_filename, copy_options);
      ++i;
      equiv_dir = fs::path(std::to_string(i));
    }

    // prototype
    std::string src_filename = clexulator_basename + "_" + bset_name + ".cc";
    if (fs::exists(src_dir / src_filename)) {
      fs::copy_file(src_dir / src_filename, dest_dir / src_filename,
                    copy_options);
    }
  }

  /// \brief Copy local clexulator and eci to test_dir and update input json
  /// with location
  ///
  /// Notes:
  /// - Assumes Clexulator files are:
  ///   - Prototype:
  ///   basis_sets/bset.<bset_name>/<project_name>_Clexulator_<bset_name>.cc
  ///   - Equivalents:
  ///   basis_sets/bset.<bset_name>/<i>/<project_name>_Clexulator_<bset_name>_<i>.cc
  ///   - Equivalents info:
  ///   basis_sets/bset.<bset_name>/equivalents_info.json
  /// - Assumes ECI files are: events/event.<bset_name>/eci.json
  /// - Assumes Event files are: events/event.<bset_name>/event.json
  void set_local_clex(std::string clex_name, std::string bset_name,
                      fs::path eci_relpath) {
    fs::path source_relpath =
        fs::path("basis_sets") / ("bset." + bset_name) /
        (project_name + "_Clexulator_" + bset_name + ".cc");
    fs::path equivalents_info_relpath = fs::path("basis_sets") /
                                        ("bset." + bset_name) /
                                        "equivalents_info.json";

    copy_local_clexulator(test_data_dir / "basis_sets", test_dir / "basis_sets",
                          bset_name, project_name + "_Clexulator");
    fs::create_directories(test_dir / eci_relpath.parent_path());
    fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                  copy_options);
    fs::copy_file(test_data_dir / equivalents_info_relpath,
                  test_dir / equivalents_info_relpath, copy_options);

    json["kwargs"]["system"]["local_clex"][clex_name]["source"] =
        (test_dir / source_relpath).string();
    json["kwargs"]["system"]["local_clex"][clex_name]["coefficients"] =
        (test_dir / eci_relpath).string();
    json["kwargs"]["system"]["local_clex"][clex_name]["equivalents_info"] =
        (test_dir / equivalents_info_relpath).string();
  }

  void set_event(std::string event_name, std::string kra_clex_name,
                 std::string freq_clex_name) {
    auto &j = json["kwargs"]["system"]["events"][event_name];

    fs::path event_relpath =
        test_dir / fs::path("events") / ("event." + event_name) / "event.json";
    fs::copy_file(test_data_dir / event_relpath, test_dir / event_relpath,
                  copy_options);

    j["event"] = event_relpath.string();
    j["kra"] = kra_clex_name;
    j["freq"] = freq_clex_name;
  }

  void write_input() { json.write(test_dir / "input.json"); }
};
