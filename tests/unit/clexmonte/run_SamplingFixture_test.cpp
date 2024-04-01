#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/metropolis.hh"
#include "casm/monte/run_management/SamplingFixture.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;
using namespace CASM;
using namespace CASM::monte;
using namespace CASM::clexmonte;

class run_SamplingFixtureTest : public test::ZrOTestSystem {
 public:
  typedef std::mt19937_64 engine_type;

  run_SamplingFixtureTest()
      : calculator(std::make_shared<canonical::Canonical_mt19937_64>(system)) {
    sampling_functions =
        canonical::Canonical_mt19937_64::standard_sampling_functions(
            calculator);
    json_sampling_functions =
        canonical::Canonical_mt19937_64::standard_json_sampling_functions(
            calculator);
    engine = std::make_shared<engine_type>();
  }

  std::shared_ptr<canonical::Canonical_mt19937_64> calculator;
  StateSamplingFunctionMap sampling_functions;
  jsonStateSamplingFunctionMap json_sampling_functions;
  ResultsAnalysisFunctionMap<config_type, statistics_type> analysis_functions;
  std::shared_ptr<engine_type> engine;
};

/// Test state sampling, using canonical Monte Carlo
TEST_F(run_SamplingFixtureTest, Test1) {
  fs::path test_dir =
      fs::current_path() / "CASM_test_projects" / "SamplingFixtureTest";
  fs::create_directories(test_dir);

  // Create conditions
  ValueMap init_conditions =
      canonical::make_conditions(600.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  state_type default_state(make_default_configuration(*system, T));
  for (Index i = 0; i < volume; ++i) {
    get_occupation(default_state)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_prim_basicstructure(*system), T};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  auto &log = CASM::log();

  // Loop over states
  for (Index i = 0; i < 8; ++i) {
    log.restart_clock();

    // Create state
    state_type state(default_state.configuration, init_conditions);
    state.conditions.scalar_values.at("temperature") = 300.0 + i * 100.0;

    std::shared_ptr<Conditions> conditions = make_conditions(*system, state);

    OccLocation occ_location(convert, occ_candidate_list);
    occ_location.initialize(get_occupation(state));
    CountType steps_per_pass = occ_location.mol_size();

    // Make potential energy calculator & set for particular supercell
    // Note: usually this happens in calculator->run
    calculator->potential =
        std::make_shared<canonical::CanonicalPotential>(system);
    calculator->potential->set(&state, conditions);
    calculator->formation_energy = calculator->potential->formation_energy();

    // Set calculator pointers
    calculator->state = &state;
    calculator->conditions = conditions;

    // Make SamplingFixture
    std::string label = "thermo";

    // Set sampling params
    SamplingParams sampling_params;
    sampling_params.sampler_names.push_back("mol_composition");
    sampling_params.sampler_names.push_back("formation_energy_corr");
    sampling_params.sampler_names.push_back("formation_energy");

    // Set analysis names
    std::vector<std::string> analysis_names = {
        "heat_capacity", "mol_susc", "param_susc", "mol_thermochem_susc",
        "param_thermochem_susc"};

    // Set completion check params
    CompletionCheckParams<statistics_type> completion_check_params;
    completion_check_params.equilibration_check_f = default_equilibration_check;
    completion_check_params.calc_statistics_f = BasicStatisticsCalculator();
    completion_check_params.cutoff_params.max_count = 100;

    // ResultsIO
    std::unique_ptr<results_io_type> results_io = nullptr;

    // Logging
    MethodLog method_log;

    SamplingFixtureParams<config_type, statistics_type> sampling_fixture_params(
        label, sampling_functions, json_sampling_functions, analysis_functions,
        sampling_params, completion_check_params, analysis_names,
        std::move(results_io), method_log);

    SamplingFixture<config_type, statistics_type, engine_type> sampling_fixture(
        sampling_fixture_params, engine);

    // Main loop
    OccEvent event;
    std::vector<Index> linear_site_index;
    std::vector<int> new_occ;
    double beta =
        1.0 / (CASM::KB * state.conditions.scalar_values.at("temperature"));
    RandomNumberGenerator<std::mt19937_64> random_number_generator;

    sampling_fixture.initialize(steps_per_pass);
    sampling_fixture.sample_data_by_count_if_due(state);
    while (!sampling_fixture.is_complete()) {
      propose_canonical_event(event, occ_location, canonical_swaps,
                              random_number_generator);

      double delta_potential_energy =
          calculator->potential->occ_delta_per_supercell(
              event.linear_site_index, event.new_occ);

      // Accept or reject event
      bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                          random_number_generator);

      // Apply accepted event
      if (accept) {
        occ_location.apply(event, get_occupation(state));
      }

      sampling_fixture.increment_step();
      sampling_fixture.sample_data_by_count_if_due(state);
    }  // main loop

    std::stringstream ss;
    ss << "samplers: " << std::endl;
    for (auto const &f : sampling_fixture.results().samplers) {
      auto const &name = f.first;
      auto const &sampler = *f.second;
      ss << name << ":" << std::endl;
      ss << "component_names: " << sampler.component_names() << std::endl;
      ss << "n_samples: " << sampler.n_samples() << std::endl;
      ss << "value: \n" << sampler.values() << std::endl;
      EXPECT_EQ(sampler.n_samples(), 100);
    }
    // std::cout << ss.str();
  }  // loop over states
}

/// Test state sampling, using canonical Monte Carlo
///
/// Test using defaults as much as possible
TEST_F(run_SamplingFixtureTest, Test2) {
  // Create conditions
  ValueMap init_conditions =
      canonical::make_conditions(600.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  state_type default_state(make_default_configuration(*system, T));
  for (Index i = 0; i < volume; ++i) {
    get_occupation(default_state)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_prim_basicstructure(*system), T};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  auto &log = CASM::log();

  // Loop over states
  for (Index i = 0; i < 8; ++i) {
    log.restart_clock();

    // Create state
    state_type state(default_state.configuration, init_conditions);
    state.conditions.scalar_values.at("temperature") = 300.0 + i * 100.0;

    std::shared_ptr<Conditions> conditions = make_conditions(*system, state);

    OccLocation occ_location(convert, occ_candidate_list);
    occ_location.initialize(get_occupation(state));
    CountType steps_per_pass = occ_location.mol_size();

    // Make potential energy calculator & set for particular supercell
    calculator->potential =
        std::make_shared<canonical::CanonicalPotential>(system);
    calculator->potential->set(&state, conditions);
    calculator->formation_energy = calculator->potential->formation_energy();

    // Set calculator pointers
    calculator->state = &state;
    calculator->conditions = conditions;

    // Make SamplingFixture
    std::string label = "thermo";

    // Set sampling params
    SamplingParams sampling_params;
    sampling_params.sampler_names.push_back("mol_composition");
    sampling_params.sampler_names.push_back("formation_energy_corr");
    sampling_params.sampler_names.push_back("formation_energy");

    // Set analysis names
    std::vector<std::string> analysis_names = {
        "heat_capacity", "mol_susc", "param_susc", "mol_thermochem_susc",
        "param_thermochem_susc"};

    // Set completion check params
    CompletionCheckParams<statistics_type> completion_check_params;
    completion_check_params.cutoff_params.max_count = 100;

    // ResultsIO
    fs::path output_dir = test_dir / "output" / label;
    bool write_trajectory = true;
    bool write_observations = true;
    std::unique_ptr<results_io_type> results_io = nullptr;

    // Logging
    MethodLog method_log;

    SamplingFixtureParams<config_type, statistics_type> sampling_fixture_params(
        label, sampling_functions, json_sampling_functions, analysis_functions,
        sampling_params, completion_check_params, analysis_names,
        std::move(results_io), method_log);

    SamplingFixture<config_type, statistics_type, engine_type> sampling_fixture(
        sampling_fixture_params, engine);

    // Main loop
    OccEvent event;
    std::vector<Index> linear_site_index;
    std::vector<int> new_occ;
    double beta =
        1.0 / (CASM::KB * state.conditions.scalar_values.at("temperature"));
    RandomNumberGenerator<std::mt19937_64> random_number_generator;

    sampling_fixture.initialize(steps_per_pass);
    sampling_fixture.sample_data_by_count_if_due(state);
    while (!sampling_fixture.is_complete()) {
      propose_canonical_event(event, occ_location, canonical_swaps,
                              random_number_generator);

      double delta_potential_energy =
          calculator->potential->occ_delta_per_supercell(
              event.linear_site_index, event.new_occ);

      // Accept or reject event
      bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                          random_number_generator);

      // Apply accepted event
      if (accept) {
        occ_location.apply(event, get_occupation(state));
      }

      sampling_fixture.increment_step();
      sampling_fixture.sample_data_by_count_if_due(state);
    }  // main loop

    std::stringstream ss;
    ss << "samplers: " << std::endl;
    for (auto const &f : sampling_fixture.results().samplers) {
      auto const &name = f.first;
      auto const &sampler = *f.second;
      ss << name << ":" << std::endl;
      ss << "component_names: " << sampler.component_names() << std::endl;
      ss << "n_samples: " << sampler.n_samples() << std::endl;
      ss << "value: \n" << sampler.values() << std::endl;
      EXPECT_EQ(sampler.n_samples(), 100);
    }
    // std::cout << ss.str();
  }  // loop over states
}
