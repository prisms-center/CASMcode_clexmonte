#include "ZrOTestSystem.hh"
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

TEST_F(CanonicalMetropolisTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  // Create conditions
  monte::VectorValueMap conditions =
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
  OccLocation occ_location(convert, occ_candidate_list);
  occ_location.initialize(get_occupation(config));
  CountType steps_per_pass = occ_location.mol_size();

  // Make supercell-specific potential energy clex calculator
  // (equal to formation energy calculator now)
  clexulator::ClusterExpansion &potential_energy_clex_calculator =
      get_formation_energy_clex(*system_data, config);
  potential_energy_clex_calculator.set(&config.dof_values);

  // Main loop
  OccEvent event;
  double beta = 1.0 / (CASM::KB * conditions.at("temperature")(0));
  double n_unitcells = get_transformation_matrix_to_super(config).determinant();
  MTRand random_number_generator;
  CountType step = 0;
  CountType pass = 0;
  while (pass < 1000) {
    // std::cout << "Propose canonical event" << std::endl;
    propose_canonical_event(event, occ_location, canonical_swaps,
                            random_number_generator);

    // std::cout << "Calculate delta_potential_energy" << std::endl;
    OccTransform const &t = event.occ_transform[0];
    Index linear_site_index = t.l;
    int new_occ = convert.occ_index(t.asym, t.to_species);
    double delta_potential_energy =
        potential_energy_clex_calculator.occ_delta_value(linear_site_index,
                                                         new_occ);

    // Accept or reject event
    // std::cout << "Check event" << std::endl;
    bool accept = metropolis_acceptance(delta_potential_energy, beta,
                                        random_number_generator);

    // Apply accepted event
    if (accept) {
      // std::cout << "Accept event" << std::endl;
      occ_location.apply(event, get_occupation(config));
      // std::cout << get_occupation(config).transpose() << std::endl;
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
