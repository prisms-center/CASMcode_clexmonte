#include "ZrOTestSystem.hh"
#include "casm/clexmonte/canonical/CanonicalPotential.hh"
#include "casm/clexmonte/canonical/conditions.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/metropolis.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class CanonicalMetropolisTest : public test::ZrOTestSystem {};

/// Test canonical monte carlo, without state sampling
TEST_F(CanonicalMetropolisTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  // Create state
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  monte::State<Configuration> state(
      make_default_configuration(*system_data, T),
      canonical::make_conditions(600.0, system_data->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}}));

  // Set initial occupation
  for (Index i = 0; i < volume; ++i) {
    get_occupation(state.configuration)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_shared_prim(*system_data),
                      get_transformation_matrix_to_super(state.configuration)};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);
  OccLocation occ_location(convert, occ_candidate_list);
  occ_location.initialize(get_occupation(state.configuration));
  CountType steps_per_pass = occ_location.mol_size();

  // Make supercell-specific potential energy clex calculator
  // (equal to formation energy calculator now)
  canonical::CanonicalPotential potential(
      get_formation_energy_clex(*system_data, state));
  set(potential, state);

  // Main loop
  OccEvent event;
  double beta =
      1.0 / (CASM::KB * state.conditions.scalar_values.at("temperature"));
  MTRand random_number_generator;
  CountType step = 0;
  CountType pass = 0;
  while (pass < 1000) {
    // std::cout << "Propose canonical event" << std::endl;
    propose_canonical_event(event, occ_location, canonical_swaps,
                            random_number_generator);

    // std::cout << "Calculate delta_potential_energy" << std::endl;
    double delta_potential_energy = potential.occ_delta_extensive_value(
        event.linear_site_index, event.new_occ);

    // Accept or reject event
    // std::cout << "Check event" << std::endl;
    bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                        random_number_generator);

    // Apply accepted event
    if (accept) {
      // std::cout << "Accept event" << std::endl;
      occ_location.apply(event, get_occupation(state.configuration));
      // std::cout << get_occupation(state.configuration).transpose() <<
      // std::endl;
    } else {
      // std::cout << "Reject event" << std::endl;
    }

    ++step;
    if (step == steps_per_pass) {
      ++pass;
      step = 0;
    }
    // std::cout << step << " " << pass << std::endl;
  }
}
