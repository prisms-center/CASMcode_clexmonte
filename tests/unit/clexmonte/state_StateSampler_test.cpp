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
  for (Index i = 0; i < 101; ++i) {
    std::cout << "Begin state: " << i << std::endl;

    // Create state
    std::cout << "Create state" << std::endl;
    State<Configuration> state(config, init_conditions);
    state.conditions.at("temperature")(0) = 300.0 + i * 10.0;

    OccLocation occ_location(convert, occ_candidate_list);
    occ_location.initialize(get_occupation(state.configuration));
    CountType steps_per_pass = occ_location.mol_size();

    // Make supercell-specific potential energy clex calculator
    // (equal to formation energy calculator now)
    std::cout << "Create potential calculator" << std::endl;
    canonical::CanonicalPotential potential(
        get_formation_energy_clex(*system_data, state));
    set(potential, state);

    // Make StateSampler
    std::vector<StateSamplingFunction<Configuration>> functions;
    functions.push_back(sampling_functions.at("comp_n"));

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
    std::cout << "Prepare for main loop" << std::endl;
    OccEvent event;
    std::vector<Index> linear_site_index;
    std::vector<int> new_occ;
    double beta = 1.0 / (CASM::KB * state.conditions.at("temperature")(0));
    MTRand random_number_generator;
    CountType step = 0;
    CountType pass = 0;
    // std::cout << "Begin main loop" << std::endl;
    // while (state_sampler.pass < 1000) {
    while (pass < 10) {
      // std::cout << "Propose canonical event" << std::endl;
      propose_canonical_event(event, occ_location, canonical_swaps,
                              random_number_generator);

      // std::cout << "Calculate delta_potential_energy" << std::endl;
      Index event_size = event.occ_transform.size();
      linear_site_index.resize(event_size);
      new_occ.resize(event_size);
      for (Index i = 0; i < event_size; ++i) {
        OccTransform const &t = event.occ_transform[i];
        linear_site_index[i] = t.l;
        new_occ[i] = convert.occ_index(t.asym, t.to_species);
        // std::cout << "linear_site_index: " << linear_site_index[i] << "/" <<
        // get_occupation(state.configuration).size() << std::endl;
        // Eigen::VectorXd const &point_corr =
        //     potential_energy_clex_calculator.correlations().point(
        //         linear_site_index[i]);
        // std::cout << "point_corr: " << point_corr.transpose() << std::endl;
      }
      // std::cout << "l: " << linear_site_index << " o: " << new_occ <<
      // std::endl; std::cout << "dof_values: " <<
      // potential_energy_clex_calculator.get() << std::endl;
      double delta_potential_energy =
          potential.occ_delta_extensive_value(linear_site_index, new_occ);

      // double delta_potential_energy =
      //     potential_energy_clex_calculator.extensive_value();
      // delta_potential_energy *= random_number_generator.rand53() - 0.5;

      // double delta_potential_energy = random_number_generator.rand53() - 0.5;

      // Accept or reject event
      // std::cout << "Check event: " << delta_potential_energy << std::endl;
      bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                          random_number_generator);

      // Apply accepted event
      if (accept) {
        // std::cout << "Accept event" << std::endl;
        occ_location.apply(event, get_occupation(state.configuration));
        // std::cout << get_occupation(config).transpose() << std::endl;
      } else {
        // std::cout << "Reject event" << std::endl;
      }

      // state_sampler.increment_step();
      // state_sampler.sample_data_if_due(state);

      ++step;
      if (step == steps_per_pass) {
        ++pass;
        step = 0;
      }
      // std::cout << std::endl;
    }  // main loop

    // std::cout << "samplers: " << std::endl;
    // for (auto const &f : state_sampler.samplers) {
    //   auto const &name = f.first;
    //   auto const &sampler = *f.second;
    //   std::cout << name << ":" << std::endl;
    //   std::cout << "component_names: " << sampler.component_names() <<
    //   std::endl; std::cout << "n_samples: " << sampler.n_samples() <<
    //   std::endl; std::cout << "value: \n" << sampler.values() << std::endl;
    // }
  }  // loop over states
}
