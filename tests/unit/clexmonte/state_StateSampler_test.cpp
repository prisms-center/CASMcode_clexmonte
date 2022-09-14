#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/canonical.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/metropolis.hh"
#include "casm/monte/state/StateSampler.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;
using namespace CASM;
using namespace CASM::monte;
using namespace CASM::clexmonte;

class state_StateSamplerTest : public test::ZrOTestSystem {
 public:
  state_StateSamplerTest()
      : calculator(
            std::make_shared<canonical::Canonical<std::mt19937_64>>(system)) {
    sampling_functions = canonical::standard_sampling_functions(calculator);
  }

  std::shared_ptr<canonical::Canonical<std::mt19937_64>> calculator;
  StateSamplingFunctionMap<Configuration> sampling_functions;
};

/// Test state sampling, using canonical Monte Carlo
TEST_F(state_StateSamplerTest, Test1) {
  // Create conditions
  monte::ValueMap init_conditions =
      canonical::make_conditions(600.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  monte::State<Configuration> default_state(
      make_default_configuration(*system, T));
  for (Index i = 0; i < volume; ++i) {
    get_occupation(default_state)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_prim_basicstructure(*system), T};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  // Loop over states
  for (Index i = 0; i < 8; ++i) {
    // Create state
    State<Configuration> state(default_state.configuration, init_conditions);
    state.conditions.scalar_values.at("temperature") = 300.0 + i * 100.0;

    OccLocation occ_location(convert, occ_candidate_list);
    occ_location.initialize(get_occupation(state));
    CountType steps_per_pass = occ_location.mol_size();

    // Make supercell-specific potential energy clex calculator
    // (equal to formation energy calculator now)
    canonical::CanonicalPotential potential(
        get_clex(*system, state, "formation_energy"));

    // Make StateSampler
    std::vector<StateSamplingFunction<Configuration>> functions;
    // for (auto f : sampling_functions) {
    //   std::cout << f.first << std::endl;
    // }
    functions.push_back(sampling_functions.at("mol_composition"));
    functions.push_back(sampling_functions.at("temperature"));
    functions.push_back(sampling_functions.at("formation_energy_corr"));
    functions.push_back(sampling_functions.at("formation_energy"));
    // no: functions.push_back(sampling_functions.at("potential_energy"));

    SAMPLE_METHOD sample_method = SAMPLE_METHOD::LINEAR;
    double sample_begin = 0.0;
    double sampling_period = 1.0;
    double samples_per_period = 1.0;
    double log_sampling_shift = 0.0;
    bool do_sample_trajectory = false;

    StateSampler state_sampler(SAMPLE_MODE::BY_PASS, functions, sample_method,
                               sample_begin, sampling_period,
                               samples_per_period, log_sampling_shift,
                               do_sample_trajectory);
    state_sampler.reset(steps_per_pass);
    state_sampler.sample_data_if_due(state);

    // Main loop
    OccEvent event;
    std::vector<Index> linear_site_index;
    std::vector<int> new_occ;
    double beta =
        1.0 / (CASM::KB * state.conditions.scalar_values.at("temperature"));
    monte::RandomNumberGenerator<std::mt19937_64> random_number_generator;

    while (state_sampler.pass < 100) {
      propose_canonical_event(event, occ_location, canonical_swaps,
                              random_number_generator);

      double delta_potential_energy = potential.occ_delta_extensive_value(
          event.linear_site_index, event.new_occ);

      // Accept or reject event
      bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                          random_number_generator);

      // Apply accepted event
      if (accept) {
        occ_location.apply(event, get_occupation(state));
      }

      state_sampler.increment_step();
      state_sampler.sample_data_if_due(state);
    }  // main loop

    std::stringstream ss;
    ss << "samplers: " << std::endl;
    for (auto const &f : state_sampler.samplers) {
      auto const &name = f.first;
      auto const &sampler = *f.second;
      ss << name << ":" << std::endl;
      ss << "component_names: " << sampler.component_names() << std::endl;
      ss << "n_samples: " << sampler.n_samples() << std::endl;
      ss << "value: \n" << sampler.values() << std::endl;
    }
    // std::cout << ss.str();
  }  // loop over states
}
