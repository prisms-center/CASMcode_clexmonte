#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/canonical/CanonicalPotential.hh"
#include "casm/clexmonte/canonical/conditions.hh"
#include "casm/clexmonte/canonical/sampling_functions.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
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

class StateSamplerTest : public test::ZrOTestSystem {
 public:
  StateSamplerTest() {
    sampling_functions = canonical::make_sampling_functions(
        system_data, canonical::canonical_tag());
  }

  StateSamplingFunctionMap<Configuration> sampling_functions;
};

/// Test state sampling, using canonical Monte Carlo
TEST_F(StateSamplerTest, Test1) {
  // Create conditions
  monte::VectorValueMap init_conditions =
      canonical::make_conditions(600.0, system_data->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  Configuration config = make_default_configuration(*system_data, T);
  for (Index i = 0; i < volume; ++i) {
    get_occupation(config)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_shared_prim(*system_data),
                      get_transformation_matrix_to_super(config)};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);

  // Loop over states
  for (Index i = 0; i < 8; ++i) {
    // Create state
    State<Configuration> state(config, init_conditions);
    state.conditions.at("temperature")(0) = 300.0 + i * 100.0;

    OccLocation occ_location(convert, occ_candidate_list);
    occ_location.initialize(get_occupation(state.configuration));
    CountType steps_per_pass = occ_location.mol_size();

    // Make supercell-specific potential energy clex calculator
    // (equal to formation energy calculator now)
    canonical::CanonicalPotential potential(
        get_formation_energy_clex(*system_data, state));
    set(potential, state);

    // Make StateSampler
    std::vector<StateSamplingFunction<Configuration>> functions;
    // for (auto f : sampling_functions) {
    //   std::cout << f.first << std::endl;
    // }
    functions.push_back(sampling_functions.at("comp_n"));
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
    double beta = 1.0 / (CASM::KB * state.conditions.at("temperature")(0));
    MTRand random_number_generator;

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
        occ_location.apply(event, get_occupation(state.configuration));
      }

      state_sampler.increment_step();
      state_sampler.sample_data_if_due(state);
    }  // main loop

    std::cout << "samplers: " << std::endl;
    for (auto const &f : state_sampler.samplers) {
      auto const &name = f.first;
      auto const &sampler = *f.second;
      std::cout << name << ":" << std::endl;
      std::cout << "component_names: " << sampler.component_names()
                << std::endl;
      std::cout << "n_samples: " << sampler.n_samples() << std::endl;
      std::cout << "value: \n" << sampler.values() << std::endl;
    }
  }  // loop over states
}
