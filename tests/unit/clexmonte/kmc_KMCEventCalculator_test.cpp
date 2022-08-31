#include "KMCCompleteEventListTestSystem.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/kmc/PrimEventCalculator.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancy_default directory
class kmc_PrimEventCalculator_Test : public KMCCompleteEventListTestSystem {};

/// \brief Test constructing event lists and calculating initial event states
///
/// Notes:
/// - FCC A-B-Va, 1NN interactions, A-Va and B-Va hops
/// - 10 x 10 x 10 (of the conventional 4-atom cell)
/// - expected runtime ~5s
TEST_F(kmc_PrimEventCalculator_Test, Test1) {
  using namespace clexmonte;
  make_prim_event_list();
  // std::cout << "#prim events: " << prim_event_list.size() << std::endl;

  // Create default state
  Eigen::Matrix3l T = test::fcc_conventional_transf_mat() * 10;
  monte::State<clexmonte::Configuration> state(
      make_default_configuration(*system, T));

  make_complete_event_list(state);
  // std::cout << "#events: " << event_list.events.size() << std::endl;

  std::vector<kmc::PrimEventCalculator> prim_event_calculators =
      kmc::make_prim_event_calculators(*system, state, prim_event_list);
  EXPECT_EQ(prim_event_calculators.size(), 24);
  // std::cout << "#prim event calculators: " << prim_event_calculators.size()
  //           << std::endl;

  // Set configuration - A, with single Va
  Eigen::VectorXi &occupation = get_occupation(state);
  occupation(0) = 2;

  // Set conditions
  monte::ValueMap conditions_map;
  conditions_map.scalar_values["temperature"] = 600.0;
  kmc::Conditions conditions(conditions_map);
  // std::cout << "beta: " << conditions.beta << std::endl;
  double expected_Ekra = 1.0;
  double expected_freq = 1e12;
  double expected_rate = expected_freq * exp(-conditions.beta * expected_Ekra);

  Index i = 0;
  Index n_allowed = 0;
  EventState event_state;
  for (auto const &event : event_list.events) {
    auto const &event_id = event.first;
    auto const &event_data = event.second;
    auto const &prim_event_data = prim_event_list[event_id.prim_event_index];
    auto const &prim_event_calculator =
        prim_event_calculators[event_id.prim_event_index];
    prim_event_calculator.calculate_event_state(event_state, conditions,
                                                event_data, prim_event_data);

    if (event_state.is_allowed) {
      // std::cout << "--- " << i << " ---" << std::endl;
      // print(std::cout, event_state, event_data, prim_event_data);
      // std::cout << std::endl;
      EXPECT_TRUE(CASM::almost_equal(event_state.dE_final, 0.0));
      EXPECT_TRUE(CASM::almost_equal(event_state.dE_activated, expected_Ekra));
      EXPECT_TRUE(CASM::almost_equal(event_state.Ekra, expected_Ekra));
      EXPECT_TRUE(CASM::almost_equal(event_state.freq, expected_freq));
      EXPECT_TRUE(CASM::almost_equal(event_state.rate, expected_rate));
      ++n_allowed;
    }
    ++i;
  }
  EXPECT_EQ(n_allowed, 12);
}
