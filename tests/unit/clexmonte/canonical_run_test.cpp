#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/io/InputData.hh"
#include "casm/clexmonte/canonical/io/json/InputData_json_io.hh"
#include "gtest/gtest.h"
#include "misc.hh"
#include "testdir.hh"

// // Test0
// #include "casm/clexmonte/system/System.hh"
// #include "casm/clexmonte/system/io/json/System_json_io.hh"
// #include "casm/clexmonte/state/sampling_functions.hh"
// #include "casm/monte/Conversions.hh"
// #include "casm/monte/checks/CompletionCheck.hh"
// #include "casm/monte/events/OccCandidate.hh"
// #include "casm/monte/methods/canonical.hh"

using namespace CASM;

// TEST(canonical_run_test, Test0) {
//   using namespace CASM::clexmonte;
//   fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
//   fs::path clexulator_src_relpath =
//       fs::path("basis_sets") / "bset.formation_energy" / "ZrO_Clexulator.cc";
//   fs::path eci_relpath = "formation_energy_eci.json";
//   fs::path output_dir_relpath = "output";
//
//   test::TmpDir test_dir;
//   test_dir.do_not_remove_on_destruction();
//   fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
//   fs::copy_file(test_data_dir / clexulator_src_relpath,
//                 test_dir / clexulator_src_relpath);
//   fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath);
//   fs::create_directories(test_dir / output_dir_relpath);
//
//   jsonParser json(test_data_dir / "input_1.json");
//   json["kwargs"]["system"]["formation_energy"]["source"] =
//       (test_dir / clexulator_src_relpath).string();
//   json["kwargs"]["system"]["formation_energy"]["coefficients"] =
//       (test_dir / eci_relpath).string();
//   json["kwargs"]["results_io"]["kwargs"]["output_dir"] =
//       (test_dir / output_dir_relpath).string();
//
//   // ---
//
//   InputParser<clexmonte::System> parser(json["kwargs"]["system"]);
//   std::runtime_error error_if_invalid{
//       "Error reading canonical Monte Carlo JSON input"};
//   report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
//   EXPECT_TRUE(parser.valid());
//   std::shared_ptr<clexmonte::System> system =
//   std::move(parser.value);
//
//   monte::VectorValueMap conditions;
//   conditions["temperature"] = Eigen::VectorXd(1);
//   conditions["temperature"](0) = 300.0;
//   conditions["mol_composition"] = Eigen::VectorXd(3);
//   conditions["mol_composition"](0) = 2.0; // Zr
//   conditions["mol_composition"](1) = 1.0; // Va
//   conditions["mol_composition"](2) = 1.0; // O
//
//   monte::State<clexmonte::Configuration> initial_state(
//       clexmonte::make_default_configuration(
//           *system, Eigen::Matrix3l::Identity()*2),
//       conditions);
//
//   // set half of sites to O
//   for (Index i=0; i<8; ++i) {
//     get_occupation(initial_state)(16+i) = 1;
//   }
//
//   // Make supercell-specific potential energy clex calculator
//   // (equal to formation energy calculator now)
//   clexulator::ClusterExpansion &potential_energy_calculator =
//       *get_clex(*system, initial_state, "formation_energy");
//
//   // Prepare supercell-specific index conversions
//   monte::Conversions convert{
//       *get_shared_prim(*system),
//       get_transformation_matrix_to_super(initial_state.configuration)};
//
//   // Prepare list of allowed swaps -- currently using all allowed
//   monte::OccCandidateList occ_candidate_list{convert};
//   std::vector<monte::OccSwap> canonical_swaps =
//       make_canonical_swaps(convert, occ_candidate_list);
//
//   MTRand random_number_generator;
//
//   std::vector<monte::StateSamplingFunction<clexmonte::Configuration>>
//   functions = {
//       make_temperature_f(system), make_comp_n_f(system),
//       make_comp_x_f(system), make_formation_energy_corr_f(system),
//       make_formation_energy_f(system),
//       make_potential_energy_f(system)};
//   monte::StateSampler<clexmonte::Configuration> state_sampler(
//       monte::SAMPLE_MODE::BY_PASS, functions);
//
//   monte::CompletionCheckParams completion_check_params;
//   completion_check_params.cutoff_params.max_count = 10;
//   monte::CompletionCheck completion_check(completion_check_params);
//
//   for (Index j=0; j<100; ++j) {
//     std::cout << "j: " << j << " temperature: " << 300 + j*10 << std::endl;
//     monte::State<clexmonte::Configuration> state = initial_state;
//     state.conditions.at("temperature")(0) = 300 + j*10;
//     // auto result = monte::canonical(initial_state,
//     potential_energy_calculator, convert, canonical_swaps,
//     random_number_generator, state_sampler, completion_check);
//     monte::OccCandidateList occ_candidate_list(convert);
//     monte::OccLocation occ_location(convert, occ_candidate_list);
//     occ_location.initialize(get_occupation(state));
//     monte::CountType steps_per_pass = occ_location.mol_size();
//
//     for (Index k=0; k<10*steps_per_pass; ++k) {
//       set(potential_energy_calculator, initial_state);
//       double epot = potential_energy_calculator.intensive_value();
//     }
//   }
//
// }

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
