#include "KMCTestSystem.hh"
#include "gtest/gtest.h"

// impact table & event lists
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/events/OccLocation.hh"

// write neighborhoods to json:
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancy_default directory
class events_impact_table_Test : public test::KMCTestSystem {
  void print_impact_info(
      std::vector<clexmonte::EventImpactInfo> const &prim_impact_info_list) {
    jsonParser json;
    json.put_array();
    for (clexmonte::EventImpactInfo const &impact : prim_impact_info_list) {
      jsonParser tjson;
      to_json(impact.phenomenal_sites, tjson["phenomenal_sites"]);
      to_json(impact.required_update_neighborhood,
              tjson["required_update_neighborhood"]);
      json.push_back(tjson);
    }
    std::cout << json << std::endl;
  }
};

/// \brief Impact neighborhood && Event lists test (FCC, 1NN interactions)
TEST_F(events_impact_table_Test, Test1) {
  setup_input_files(false /*use_sparse_format_eci*/);

  std::vector<clexmonte::PrimEventData> prim_event_list =
      make_prim_event_list(*system);
  EXPECT_EQ(prim_event_list.size(), 24);

  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list =
      make_prim_impact_info_list(*system, prim_event_list,
                                 {"formation_energy"});

  EXPECT_EQ(prim_impact_info_list.size(), 24);
  for (auto const &impact : prim_impact_info_list) {
    EXPECT_EQ(impact.required_update_neighborhood.size(), 20);
  }

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  monte::State<clexmonte::Configuration> state(
      make_default_configuration(*system, T));
  monte::OccLocation occ_location{get_index_conversions(*system, state),
                                  get_occ_candidate_list(*system, state)};
  occ_location.initialize(get_occupation(state));

  clexmonte::CompleteEventList event_list = clexmonte::make_complete_event_list(
      prim_event_list, prim_impact_info_list, occ_location);
  EXPECT_EQ(event_list.impact_table.size(), 1000 * 24);

  for (auto const &impacted : event_list.impact_table) {
    EXPECT_EQ(impacted.second.size(), 708);
  }
  EXPECT_EQ(event_list.events.size(), 1000 * 24);
}

/// \brief Simpler test:
/// - Single A-Va event,
/// - Only point formation energy contribution
/// - Only constant KRA,freq contribution
/// --> 45 impacted events:
///   - 10 * 4 (unaligned events)
///     + 1 * 3 (aligned reverse events)
///     + 1 * 2 (aligned equivalent events)
TEST_F(events_impact_table_Test, Test2) {
  setup_input_files(false /*use_sparse_format_eci*/);

  clexulator::SparseCoefficients constant_eci;
  constant_eci.index = {0};
  constant_eci.value = {0.};

  clexulator::SparseCoefficients point_eci;
  point_eci.index = {0, 1};
  point_eci.value = {0., 0.};

  system->clex_data.at("formation_energy").coefficients = point_eci;
  system->local_multiclex_data.at("A_Va_1NN").coefficients[0] = constant_eci;
  system->local_multiclex_data.at("A_Va_1NN").coefficients[1] = constant_eci;
  system->event_type_data.erase("B_Va_1NN");

  // 6 NN clusters, forward & reverse -> 12
  std::vector<clexmonte::PrimEventData> prim_event_list =
      make_prim_event_list(*system);
  EXPECT_EQ(prim_event_list.size(), 12);

  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list =
      make_prim_impact_info_list(*system, prim_event_list,
                                 {"formation_energy"});
  EXPECT_EQ(prim_impact_info_list.size(), 12);

  // only the event sites themselves -> 2
  for (auto const &impact : prim_impact_info_list) {
    EXPECT_EQ(impact.required_update_neighborhood.size(), 2);
  }

  // print_impact_info(prim_impact_info_list);

  // Create config
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 10;
  monte::State<clexmonte::Configuration> state(
      make_default_configuration(*system, T));
  monte::OccLocation occ_location{get_index_conversions(*system, state),
                                  get_occ_candidate_list(*system, state)};
  occ_location.initialize(get_occupation(state));

  clexmonte::CompleteEventList event_list = clexmonte::make_complete_event_list(
      prim_event_list, prim_impact_info_list, occ_location);
  EXPECT_EQ(event_list.impact_table.size(), 1000 * 12);

  for (auto const &impacted : event_list.impact_table) {
    EXPECT_EQ(impacted.second.size(), 46);
  }
  EXPECT_EQ(event_list.events.size(), 1000 * 12);
}

// /// \brief Useful for big supercell tests
// TEST_F(events_impact_table_Test, Test3) {
//
//   std::cout << "Construct prim event list... " << std::endl;
//   std::vector<clexmonte::PrimEventData> prim_event_list =
//       make_prim_event_list(*system);
//   std::cout << "  Done" << std::endl;
//   EXPECT_EQ(prim_event_list.size(), 24);
//
//   std::cout << "Construct prim impact info list... " << std::endl;
//   std::vector<clexmonte::EventImpactInfo> prim_impact_info_list =
//       make_prim_impact_info_list(*system, prim_event_list,
//                                  {"formation_energy"});
//   std::cout << "  Done" << std::endl;
//   EXPECT_EQ(prim_impact_info_list.size(), 24);
//
//   for (auto const &impact : prim_impact_info_list) {
//     EXPECT_EQ(impact.required_update_neighborhood.size(), 20);
//   }
//
//   //print_impact_info(prim_impact_info_list);
//
//   // Create config
//   std::cout << "Construct supercell and state... " << std::endl;
//   Index dim = 30;
//   Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * dim;
//   monte::State<clexmonte::Configuration> state(
//       make_default_configuration(*system, T));
//   std::cout << "  Done" << std::endl;
//
//   std::cout << "Construct atom tracking lists... " << std::endl;
//   monte::OccLocation occ_location{get_index_conversions(*system, state),
//                                   get_occ_candidate_list(*system, state)};
//   occ_location.initialize(get_occupation(state));
//   std::cout << "  Done" << std::endl;
//
//   Index vol = dim * dim * dim;
//   std::cout << "Construct event list and impact table... " << std::endl;
//   clexmonte::CompleteEventList event_list =
//   clexmonte::make_complete_event_list(prim_event_list,
//   prim_impact_info_list, occ_location); std::cout << "  Done" << std::endl;
//   EXPECT_EQ(event_list.impact_table.size(), vol * 24);
//
//   for (auto const &impacted : event_list.impact_table) {
//     EXPECT_EQ(impacted.second.size(), 707);
//   }
//   EXPECT_EQ(event_list.events.size(), vol * 24);
//
//   // double sum = 0.0;
//   // while (true) {
//   //   sum += 1.0;
//   // }
// }
