#include "ZrOTestSystem.hh"
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
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class canonical_MetropolisTest : public test::ZrOTestSystem {};

/// Test canonical monte carlo, without state sampling
TEST_F(canonical_MetropolisTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  // Create state
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  Index volume = T.determinant();
  state_type state(
      make_default_configuration(*system, T),
      canonical::make_conditions(600.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}}));
  std::shared_ptr<Conditions> conditions = make_conditions(*system, state);
  // Set initial occupation
  for (Index i = 0; i < volume; ++i) {
    get_occupation(state)(2 * volume + i) = 1;
  }

  // Prepare supercell-specific index conversions
  Conversions convert{*get_prim_basicstructure(*system),
                      get_transformation_matrix_to_super(state)};
  OccCandidateList occ_candidate_list(convert);
  std::vector<OccSwap> canonical_swaps =
      make_canonical_swaps(convert, occ_candidate_list);
  OccLocation occ_location(convert, occ_candidate_list);
  occ_location.initialize(get_occupation(state));
  CountType steps_per_pass = occ_location.mol_size();

  // Make potential energy calculator & set for particular supercell
  canonical::CanonicalPotential potential(system);
  potential.set(&state, conditions);

  // Main loop
  OccEvent event;
  double beta =
      1.0 / (CASM::KB * state.conditions.scalar_values.at("temperature"));
  monte::RandomNumberGenerator<std::mt19937_64> random_number_generator;
  CountType step = 0;
  CountType pass = 0;
  std::cout << "Run";
  while (pass < 100) {
    // std::cout << "Propose canonical event" << std::endl;
    propose_canonical_event(event, occ_location, canonical_swaps,
                            random_number_generator);

    // std::cout << "Calculate delta_potential_energy" << std::endl;
    double delta_potential_energy = potential.occ_delta_per_supercell(
        event.linear_site_index, event.new_occ);

    // Accept or reject event
    // std::cout << "Check event" << std::endl;
    bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                        random_number_generator);

    // Apply accepted event
    if (accept) {
      // std::cout << "Accept event" << std::endl;
      occ_location.apply(event, get_occupation(state));
      // std::cout << get_occupation(state).transpose() <<
      // std::endl;
    } else {
      // std::cout << "Reject event" << std::endl;
    }

    ++step;
    if (step == steps_per_pass) {
      ++pass;
      std::cout << ".";
      step = 0;
    }
    // std::cout << step << " " << pass << std::endl;
  }
  std::cout << std::endl << "Done" << std::endl;
}
