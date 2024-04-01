#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/run/FixedConfigGenerator.hh"
#include "casm/clexmonte/run/IncrementalConditionsStateGenerator.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/io/json/SparseCoefficients_json_io.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/checks/CompletionCheck.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/SamplingFixture.hh"
#include "casm/monte/run_management/io/json/jsonResultsIO_impl.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"
#include "casm/monte/sampling/SamplingParams.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "gtest/gtest.h"
#include "misc.hh"
#include "testdir.hh"

using namespace CASM;

TEST(canonical_fullrun_test, Test1) {
  // Copy test input data to a temperorary directory
  fs::path test_data_dir = test::data_dir("clexmonte") / "Clex_ZrO_Occ";
  fs::path clexulator_src_relpath = fs::path("basis_sets") /
                                    "bset.formation_energy" /
                                    "ZrO_Clexulator_formation_energy.cc";
  fs::path eci_relpath = "formation_energy_eci.json";
  fs::path prim_relpath = "prim.json";

  fs::path test_dir =
      fs::current_path() / "CASM_test_projects" / "canonical_fullrun_test";
  fs::copy_options copy_options = fs::copy_options::skip_existing;
  fs::create_directories(test_dir / clexulator_src_relpath.parent_path());
  fs::copy_file(test_data_dir / clexulator_src_relpath,
                test_dir / clexulator_src_relpath, copy_options);
  fs::copy_file(test_data_dir / eci_relpath, test_dir / eci_relpath,
                copy_options);
  fs::copy_file(test_data_dir / prim_relpath, test_dir / prim_relpath,
                copy_options);

  // Set Clexulator compilation options
  //   ex: g++ -O3 -Wall -fPIC --std=c++17 -I/path/to/include
  std::string default_clexulator_compile_options =
      //
      // uses $CASM_CXX, else default="g++"
      RuntimeLibrary::default_cxx().first + " " +
      //
      // uses $CASM_CXXFLAGS, else default="-O3 -Wall -fPIC --std=c++17"
      RuntimeLibrary::default_cxxflags().first + " " +
      //
      // uses -I$CASM_INCLUDEDIR,
      //   else -I$CASM_PREFIX/include,
      //   else tries to find "ccasm" or "casm" executable on PATH and looks
      //     for standard include paths relative from there,
      //   else fails with "/not/found"
      include_path(RuntimeLibrary::default_casm_includedir().first);

  // Set Clexulator shared object compilation options
  //   ex: g++ -shared -L/path/to/lib -lcasm_global -lcasm_crystallography
  //     -lcasm_clexulator -lcasm_monte
  std::string default_clexulator_so_options =
      //
      // uses $CASM_CXX, else default="g++"
      RuntimeLibrary::default_cxx().first + " " +
      //
      // uses $CASM_SOFLAGS, else default="-shared"
      RuntimeLibrary::default_soflags().first + " " +
      //
      // uses -L$CASM_LIBDIR,
      //   else -L$CASM_PREFIX/lib,
      //   else tries to find "ccasm" or "casm" executables on PATH and looks
      //     for libcasm at standard relative paths from there,
      //   else fails with "-L/not/found"
      link_path(RuntimeLibrary::default_casm_libdir().first) + " " +
      //
      // requires libcasm_clexulator:
      "-lcasm_clexulator ";

  // Create an output directory
  fs::path output_dir_relpath = "output";
  fs::path output_dir = test_dir / output_dir_relpath;

  // Error message
  std::runtime_error error_if_invalid{
      "Error reading canonical Monte Carlo JSON input"};

  // ### Construct system data

  // - Construct prim
  jsonParser prim_json(test_dir / prim_relpath);
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      std::make_shared<xtal::BasicStructure const>(read_prim(prim_json, TOL));

  // - Construct composition::CompositionConverter
  std::vector<std::string> components = {"Zr", "Va", "O"};

  Eigen::VectorXd origin;
  origin.resize(3);
  origin << 2.0, 2.0, 0.0;

  Eigen::MatrixXd end_members;
  end_members.resize(3, 1);
  end_members.col(0) << 2.0, 0.0, 2.0;

  composition::CompositionConverter composition_converter(components, origin,
                                                          end_members);

  // - Construct system data
  std::shared_ptr<clexmonte::System> system =
      std::make_shared<clexmonte::System>(
          shared_prim,  // std::shared_ptr<xtal::BasicStructure const> const &
          composition_converter  // composition::CompositionConverter const &
      );

  // - Construct clexulator::Clexulator for formation energy
  fs::path clexulator_src = test_dir / clexulator_src_relpath;
  std::string clexulator_name = clexulator_src.stem();
  fs::path clexulator_dirpath = clexulator_src.parent_path();
  std::string clexulator_compile_options = default_clexulator_compile_options;
  std::string clexulator_so_options = default_clexulator_so_options;
  std::shared_ptr<clexulator::Clexulator> clexulator =
      std::make_shared<clexulator::Clexulator>(clexulator::make_clexulator(
          clexulator_name, clexulator_dirpath, system->prim_neighbor_list,
          clexulator_compile_options, clexulator_so_options));

  // - Construct clexulator::SparseCoefficients for formation energy
  jsonParser eci_json(test_dir / eci_relpath);
  InputParser<clexulator::SparseCoefficients> eci_parser(eci_json);
  report_and_throw_if_invalid(eci_parser, CASM::log(), error_if_invalid);
  clexulator::SparseCoefficients eci = *eci_parser.value;

  // - Add formation energy basis set to `system`
  system->basis_sets.emplace("formation_energy", clexulator);

  // - Add formation energy clex to `system`
  clexmonte::ClexData formation_energy_clex_data;
  formation_energy_clex_data.basis_set_name = "formation_energy";
  formation_energy_clex_data.coefficients = eci;
  system->clex_data.emplace("formation_energy", formation_energy_clex_data);

  // ### Construct the canonical calculator
  typedef clexmonte::canonical::Canonical_mt19937_64 calculation_type;
  typedef calculation_type::engine_type engine_type;
  auto calculation = std::make_shared<calculation_type>(system);

  // Make random number engine
  std::shared_ptr<engine_type> engine = std::make_shared<engine_type>();

  // ### Construct sampling functions & analysis functions
  monte::StateSamplingFunctionMap sampling_functions =
      calculation_type::standard_sampling_functions(calculation);
  EXPECT_EQ(sampling_functions.size(), 6);
  EXPECT_EQ(sampling_functions.count("temperature"), 1);
  EXPECT_EQ(sampling_functions.count("mol_composition"), 1);
  EXPECT_EQ(sampling_functions.count("param_composition"), 1);
  EXPECT_EQ(sampling_functions.count("formation_energy_corr"), 1);
  EXPECT_EQ(sampling_functions.count("formation_energy"), 1);
  EXPECT_EQ(sampling_functions.count("potential_energy"), 1);
  monte::jsonStateSamplingFunctionMap json_sampling_functions =
      calculation_type::standard_json_sampling_functions(calculation);
  EXPECT_EQ(json_sampling_functions.size(), 1);
  EXPECT_EQ(json_sampling_functions.count("config"), 1);
  std::map<std::string, clexmonte::results_analysis_function_type>
      analysis_functions =
          calculation_type::standard_analysis_functions(calculation);
  // - Add custom sampling functions if desired...
  // state_sampling_function_type f {
  //     "potential_energy", // sampler name
  //     "Potential energy of the state (normalized per primitive cell)", //
  //     description 1,  // number of components in "potential_energy"
  //     [system](state_type const &state) {
  //       return state.properties.at("potential_energy");
  //     });
  // sampling_functions.emplace(f.name, f);

  // ### Construct the state generator

  // - Specify the supercell transformation_matrix_to_super
  Eigen::Matrix3l transformation_matrix_to_super;
  transformation_matrix_to_super.col(0) << 10, 0, 0;
  transformation_matrix_to_super.col(1) << 0, 10, 0;
  transformation_matrix_to_super.col(2) << 0, 0, 10;

  // - Construct an initial configuration (use default DoF values)
  clexmonte::Configuration initial_configuration =
      clexmonte::make_default_configuration(*system,
                                            transformation_matrix_to_super);
  // - Construct a configuration generator
  auto config_generator = notstd::make_unique<clexmonte::FixedConfigGenerator>(
      initial_configuration);

  // - Construct initial conditions
  monte::ValueMap initial_conditions = clexmonte::canonical::make_conditions(
      300.0,                  // temperature (K)
      composition_converter,  // composition converter
      {{"Zr", 2.},            // composition values (#/unit cell)
       {"O", 2. / 6.},
       {"Va", 10. / 6.}});

  // - Construct conditions increment
  monte::ValueMap conditions_increment =
      clexmonte::canonical::make_conditions_increment(
          10.0,                   // temperature (K)
          composition_converter,  // composition converter
          {{"Zr", 0.0},           // composition values (#/unit cell)
           {"O", 0.01},
           {"Va", -0.01}});

  // - Specify number of states (includes initial conditions)
  Index n_states = 11;

  // - Specify if dependent runs
  //   (if true, use final configuration at previous state as the
  //   initial configuration for the next state)
  bool dependent_runs = true;

  // - Specify if any conditions should be treated as "dependent"
  //   - For example, instead of setting composition as a independent
  //     condition, "mol_composition" could be a calculated from
  //     the generated configuration.
  std::vector<clexmonte::StateModifyingFunction> modifiers;

  // - Construct the state generator
  clexmonte::RunDataOutputParams output_params;
  output_params.do_save_all_final_states = true;
  output_params.output_dir = output_dir;
  clexmonte::IncrementalConditionsStateGenerator state_generator(
      system, output_params, std::move(config_generator), initial_conditions,
      conditions_increment, n_states, dependent_runs, modifiers);

  // ### Construct monte::SamplingParams
  monte::SamplingParams sampling_params;

  // - Sample by step, pass, or time
  sampling_params.sample_mode = monte::SAMPLE_MODE::BY_PASS;

  // - Sample linearly or logarithmically
  //
  // Default=SAMPLE_METHOD::LINEAR
  //
  // For SAMPLE_METHOD::LINEAR, take the n-th (n=0,1,2,...) sample when:
  //
  //    sample/pass = round( begin + period * n )
  //           time = begin + period * n
  //
  // For SAMPLE_METHOD::LOG, take the n-th sample when:
  //
  //    sample/pass = round( begin + base ^ (n + shift)
  //           time = begin + base ^ (n + shift)
  //
  sampling_params.sample_method = monte::SAMPLE_METHOD::LINEAR;
  sampling_params.begin = 1.0;
  sampling_params.period = 1.0;
  sampling_params.base = 1.0;
  sampling_params.shift = 0.0;

  // - What sampling functions to sample
  sampling_params.sampler_names = std::vector<std::string>(
      {"temperature", "mol_composition", "param_composition",
       "formation_energy_corr", "formation_energy", "potential_energy"});

  sampling_params.json_sampler_names = std::vector<std::string>({"config"});

  std::vector<std::string> analysis_names = {
      "heat_capacity", "mol_susc", "param_susc", "mol_thermochem_susc",
      "param_thermochem_susc"};

  // - Store configurations at sampling time
  sampling_params.do_sample_trajectory = false;

  // ### Construct monte::CompletionCheckParams
  monte::CompletionCheckParams<clexmonte::statistics_type>
      completion_check_params;
  completion_check_params.equilibration_check_f =
      monte::default_equilibration_check;
  completion_check_params.calc_statistics_f =
      monte::BasicStatisticsCalculator();

  // - Set monte::CutoffCheckParams
  completion_check_params.cutoff_params.min_count = std::nullopt;
  completion_check_params.cutoff_params.max_count = 100;
  completion_check_params.cutoff_params.min_sample = std::nullopt;
  completion_check_params.cutoff_params.max_sample = std::nullopt;

  // - Set requested precision for convergence
  converge(sampling_functions, completion_check_params)
      .set_abs_precision("formation_energy", 0.001)
      .set_abs_precision("formation_energy_corr", 0.01)
      .set_abs_precision("mol_composition", 0.01, {"O"});

  // - Set other completion check parameters or use defaults
  // completion_check_params.confidence = 0.95; // default=0.95
  // completion_check_params.log_spacing = true;  // default=true
  // completion_check_params.check_begin = 0.0; // default=0
  // completion_check_params.check_period = 10.0;  // default=10
  // completion_check_params.check_per_period = 1.0;  // default=1
  // completion_check_params.check_shift = 1.0;  // default=1

  // ### Construct monte::jsonResultsIO
  fs::path output_thermo_dir = output_dir / "thermo";
  bool write_trajectory = true;
  bool write_observations = true;
  auto results_io =
      std::make_unique<monte::jsonResultsIO<clexmonte::results_type>>(
          output_thermo_dir, write_trajectory, write_observations);

  // ~~~~ Run ~~~~

  // Create monte::MethodLog
  monte::MethodLog method_log;
  method_log.logfile_path = output_thermo_dir / "status.json";
  method_log.log_frequency = 60;  // seconds

  bool global_cutoff = true;

  std::vector<clexmonte::sampling_fixture_params_type> sampling_fixture_params;
  std::string label = "thermo";
  sampling_fixture_params.emplace_back(
      "thermo", sampling_functions, json_sampling_functions, analysis_functions,
      sampling_params, completion_check_params, analysis_names,
      std::move(results_io), method_log);

  clexmonte::run_series(*calculation, engine, state_generator,
                        sampling_fixture_params, global_cutoff);

  // check output/ files presence
  EXPECT_TRUE(fs::exists(output_dir / "completed_runs.json"));

  // check reading output/completed_runs.json
  clexmonte::IncrementalConditionsStateGenerator state_generator_2(
      system, output_params, std::move(config_generator), initial_conditions,
      conditions_increment, n_states, dependent_runs, modifiers);
  state_generator_2.read_completed_runs();
  EXPECT_EQ(state_generator_2.n_completed_runs(), 11);

  // check output/thermo files presence
  EXPECT_TRUE(fs::exists(output_thermo_dir));
  EXPECT_TRUE(fs::exists(output_thermo_dir / "summary.json"));
  EXPECT_EQ(test::file_count(output_thermo_dir), 13);
  for (int i = 1; i <= 11; ++i) {
    fs::path run_dir =
        output_thermo_dir / (std::string("run.") + std::to_string(i));
    EXPECT_TRUE(fs::exists(run_dir));
    EXPECT_EQ(test::file_count(run_dir), 2);
    EXPECT_TRUE(fs::exists(run_dir / "observations.json"));
    EXPECT_TRUE(fs::exists(run_dir / "trajectory.json"));
  }

  // remove output/thermo/ files
  if (fs::exists(output_thermo_dir / "status.json")) {
    fs::remove(output_thermo_dir / "status.json");
  }
  for (int i = 1; i <= 11; ++i) {
    fs::path run_dir =
        output_thermo_dir / (std::string("run.") + std::to_string(i));
    fs::remove(run_dir / "observations.json");
    fs::remove(run_dir / "trajectory.json");
    fs::remove(run_dir);
  }
  fs::remove(output_thermo_dir / "summary.json");
  fs::remove(output_thermo_dir);

  // remove output/ files
  fs::remove(output_dir / "completed_runs.json");

  // remote output/
  fs::remove(output_dir);
}
