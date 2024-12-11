#include "KMCTestSystem.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

namespace test {

/// \brief Default test project - if you use this, don't copy other
///     files into the test project, just use the default test fixture
KMCTestSystem::KMCTestSystem()
    : KMCTestSystem(
          "FCC_binary_vacancy", "FCCBinaryVacancy_default",
          test::data_dir("clexmonte") / "kmc" / "system_template.json") {}

void KMCTestSystem::setup_input_files(bool use_sparse_format_eci) {
  if (!use_sparse_format_eci) {
    set_clex("formation_energy", "default", "formation_energy_eci.json");

    {
      fs::path event_relpath = fs::path("kmc_events") / "event.A_Va_1NN";
      set_local_basis_set("A_Va_1NN");
      set_event("A_Va_1NN", event_relpath / "kra_eci.json",
                event_relpath / "freq_eci.json");
    }

    {
      fs::path event_relpath = fs::path("kmc_events") / "event.B_Va_1NN";
      set_local_basis_set("B_Va_1NN");
      set_event("B_Va_1NN", event_relpath / "kra_eci.json",
                event_relpath / "freq_eci.json");
    }
  } else {
    set_clex("formation_energy", "default", "formation_energy_sparse_eci.json");

    {
      fs::path event_relpath = fs::path("kmc_events") / "event.A_Va_1NN";
      set_local_basis_set("A_Va_1NN");
      set_event("A_Va_1NN", event_relpath / "kra_sparse_eci.json",
                event_relpath / "freq_sparse_eci.json");
    }

    {
      fs::path event_relpath = fs::path("kmc_events") / "event.B_Va_1NN";
      set_local_basis_set("B_Va_1NN");
      set_event("B_Va_1NN", event_relpath / "kra_sparse_eci.json",
                event_relpath / "freq_sparse_eci.json");
    }
  }
  write_input();
  make_system();
}

/// \param _project_name Name of project test files to use
///     (clexmonte/data/<_project_name>)
/// \param _test_dir_name Name of directory where test files should be
///     copied and tested (CASM_test_projects/<_test_dir_name>)
/// \param _input_file_path Path to template input file that will be
///     updated and parsed to construct the System. For example,
///     use `test::data_dir("clexmonte") / "kmc" / "system_template.json".
///     It is expected to have at least:
///     \code
///     {
///       "kwargs": {
///          "system": {
///            "prim": {...},
///            "composition_axes": {...}
///          }
///     }
///     \endcode
///     .
///
KMCTestSystem::KMCTestSystem(std::string _project_name,
                             std::string _test_dir_name,
                             fs::path _input_file_path)
    : project_name(_project_name),
      test_data_dir(test::data_dir("clexmonte") / _project_name),
      test_dir(fs::current_path() / "CASM_test_projects" / _test_dir_name),
      copy_options(fs::copy_options::skip_existing),
      json(_input_file_path) {}

/// \brief Copy formation_energy Clexulator and ECI to test_dir and update
///     input json with location
///
/// Notes:
/// - Assumes Clexulator file is:
///   - basis_sets/bset.<bset_name>/<project_name>_Clexulator_<bset_name>.cc
/// - ECI files is <eci_relpath>
void KMCTestSystem::set_clex(std::string clex_name, std::string bset_name,
                             fs::path eci_relpath) {
  fs::path clexulator_src_relpath =
      fs::path("basis_sets") / ("bset." + bset_name) /
      (project_name + "_Clexulator_" + bset_name + ".cc");
  fs::path basis_relpath =
      fs::path("basis_sets") / ("bset." + bset_name) / "basis.json";

  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::copy_file(test_data_dir / basis_relpath, test_dir / basis_relpath,
                copy_options);
  fs::create_directories(test_dir / eci_relpath.parent_path());
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);

  json["kwargs"]["system"]["basis_sets"][bset_name]["source"] =
      (test_dir / clexulator_src_relpath).string();
  json["kwargs"]["system"]["basis_sets"][bset_name]["basis"] =
      (test_dir / basis_relpath).string();
  json["kwargs"]["system"]["clex"][clex_name]["basis_set"] = bset_name;
  json["kwargs"]["system"]["clex"][clex_name]["coefficients"] =
      (test_dir / eci_relpath).string();
}

void KMCTestSystem::copy_local_clexulator(fs::path src_basis_sets_dir,
                                          fs::path dest_basis_sets_dir,
                                          std::string bset_name,
                                          std::string clexulator_basename) {
  fs::path src_dir = src_basis_sets_dir / (std::string("bset.") + bset_name);
  fs::path dest_dir = dest_basis_sets_dir / (std::string("bset.") + bset_name);

  // equivalents
  Index i = 0;
  fs::path equiv_dir = fs::path(std::to_string(i));
  while (fs::exists(src_dir / equiv_dir)) {
    std::string src_filename =
        clexulator_basename + "_" + bset_name + "_" + std::to_string(i) + ".cc";
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
///   - Basis info:
///   basis_sets/bset.<bset_name>/basis.json
/// - Assumes ECI files are: events/event.<bset_name>/eci.json
/// - Assumes Event files are: events/event.<bset_name>/event.json
void KMCTestSystem::set_local_basis_set(std::string bset_name) {
  fs::path source_relpath = fs::path("basis_sets") / ("bset." + bset_name) /
                            (project_name + "_Clexulator_" + bset_name + ".cc");
  fs::path basis_relpath =
      fs::path("basis_sets") / ("bset." + bset_name) / "basis.json";
  fs::path equivalents_info_relpath =
      fs::path("basis_sets") / ("bset." + bset_name) / "equivalents_info.json";

  copy_local_clexulator(test_data_dir / "basis_sets", test_dir / "basis_sets",
                        bset_name, project_name + "_Clexulator");
  fs::copy_file(test_data_dir / basis_relpath, test_dir / basis_relpath,
                copy_options);
  fs::copy_file(test_data_dir / equivalents_info_relpath,
                test_dir / equivalents_info_relpath, copy_options);

  json["kwargs"]["system"]["local_basis_sets"][bset_name]["basis"] =
      (test_dir / basis_relpath).string();
  json["kwargs"]["system"]["local_basis_sets"][bset_name]["source"] =
      (test_dir / source_relpath).string();
  json["kwargs"]["system"]["local_basis_sets"][bset_name]["equivalents_info"] =
      (test_dir / equivalents_info_relpath).string();
}

void KMCTestSystem::set_event(std::string event_name,
                              std::string kra_eci_relpath,
                              std::string freq_eci_relpath) {
  fs::create_directories(test_dir / "kmc_events");
  fs::copy_file(test_data_dir / "kmc_events" / "event_system.json",
                test_dir / "kmc_events" / "event_system.json", copy_options);
  json["kwargs"]["system"]["event_system"] =
      (test_dir / fs::path("kmc_events") / "event_system.json").string();

  fs::path event_relpath =
      fs::path("kmc_events") / ("event." + event_name) / "event.json";
  fs::create_directories((test_dir / event_relpath).parent_path());
  fs::copy_file(test_data_dir / event_relpath, test_dir / event_relpath,
                copy_options);

  fs::create_directories((test_dir / kra_eci_relpath).parent_path());
  fs::copy_file(test_data_dir / kra_eci_relpath, test_dir / kra_eci_relpath,
                copy_options);

  fs::create_directories((test_dir / freq_eci_relpath).parent_path());
  fs::copy_file(test_data_dir / freq_eci_relpath, test_dir / freq_eci_relpath,
                copy_options);

  auto &j = json["kwargs"]["system"]["kmc_events"][event_name];
  j["event"] = (test_dir / event_relpath).string();
  j["local_basis_set"] = event_name;
  j["coefficients"]["kra"] = (test_dir / kra_eci_relpath).string();
  j["coefficients"]["freq"] = (test_dir / freq_eci_relpath).string();
}

void KMCTestSystem::write_input() { json.write(test_dir / "input.json"); }

void KMCTestSystem::make_system() {
  EXPECT_TRUE(json.find_at(fs::path("kwargs") / "system") != json.end())
      << "Bad KMCTestSystem JSON input";
  std::vector<fs::path> search_path;
  InputParser<clexmonte::System> parser(json["kwargs"]["system"], search_path);
  std::runtime_error error_if_invalid{"Error reading KMC System JSON input"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  EXPECT_TRUE(parser.value != nullptr) << "Bad KMCTestSystem parsing";
  system = std::shared_ptr<clexmonte::System>(parser.value.release());
}

}  // namespace test
