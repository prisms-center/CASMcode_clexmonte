#include "KMCCompleteEventCalculatorTestSystem.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/kinetic/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancy_default directory
class events_RejectionFree_Test
    : public test::KMCCompleteEventCalculatorTestSystem {};

/// \brief Test constructing event lists and calculating initial event states
///
/// Notes:
/// - FCC A-B-Va, 1NN interactions, A-Va and B-Va hops
/// - 10 x 10 x 10 (of the conventional 4-atom cell)
/// - expected runtime ~8s
TEST_F(events_RejectionFree_Test, Test1) {
  using namespace clexmonte;
  setup_input_files(false /*use_sparse_format_eci*/);

  // --- State setup ---

  // Create default state
  Index dim = 5;
  Eigen::Matrix3l T = test::fcc_conventional_transf_mat() * dim;
  monte::State<clexmonte::Configuration> state(
      make_default_configuration(*system, T));

  // Set configuration - A, with single Va
  Eigen::VectorXi &occupation = get_occupation(state);
  occupation(0) = 2;
  occupation(1) = 2;
  occupation(2) = 2;

  // Set conditions
  state.conditions.scalar_values.emplace("temperature", 600);

  // --- KMC implementation ---

  // Make calculator
  make_complete_event_calculator(state);

  // Make selector
  lotto::RejectionFreeEventSelector selector(
      event_calculator,
      clexmonte::make_complete_event_id_list(T.determinant(), prim_event_list),
      event_list.impact_table);

  // Run
  std::cout << "Begin run: " << std::endl;
  double time = 0.0;
  Index n_steps = 1000;
  auto unitcell_index_converter =
      occ_location->convert().unitcell_index_converter();
  auto index_converter = occ_location->convert().index_converter();
  auto const &basicstructure = *system->prim->basicstructure;
  xtal::UnitCell translation;
  clexmonte::EventID id;
  double time_step;

  Index i = 0;

  auto print_state = [&]() {
    std::cout << "step: " << i << " simulated_time: " << time << std::endl;
    std::cout << "occupation: " << std::endl;
    for (Index l = 0; l < occupation.size(); ++l) {
      if (occupation(l) != 0) {
        std::cout << "l: " << l << "   occ: " << occupation(l)
                  << "  site: " << index_converter(l) << "   cart: "
                  << index_converter(l)
                         .coordinate(basicstructure)
                         .const_cart()
                         .transpose()
                  << std::endl;
      }
    }
  };

  for (; i < n_steps; ++i) {
    if (i % 100 == 0) {
      print_state();
    }

    std::tie(id, time_step) = selector.select_event();

    // std::cout << "--- " << i << " ---" << std::endl;
    // std::cout << "prim_event_index: " << id.prim_event_index << std::endl;
    // std::cout << "unitcell_index: " << id.unitcell_index
    //     << "   translation: " <<
    //     unitcell_index_converter(id.unitcell_index).transpose()
    //     << std::endl;
    // std::cout << std::endl;

    // Apply accepted event
    auto const &event_data = event_list.events.at(id);
    occ_location->apply(event_data.event, occupation);
    time += time_step;
  }
  print_state();
  std::cout << "Done" << std::endl << std::endl;

  EXPECT_TRUE(true);
}
