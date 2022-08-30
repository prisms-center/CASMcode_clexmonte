#include "KMCTestSystem.hh"
#include "gtest/gtest.h"

// impact table
#include "casm/clexmonte/events/event_methods.hh"

// write neighborhoods to json:
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancy_default directory
class kmc_impact_table_Test : public KMCTestSystem {};

TEST_F(kmc_impact_table_Test, Test1) {
  EXPECT_TRUE(system != nullptr);
  EXPECT_EQ(system->basis_sets.size(), 1);
  EXPECT_EQ(system->local_basis_sets.size(), 2);
  EXPECT_EQ(system->clex_data.size(), 1);
  EXPECT_EQ(system->local_multiclex_data.size(), 2);
  EXPECT_EQ(system->event_type_data.size(), 2);

  clexulator::SparseCoefficients constant_eci;
  constant_eci.index.push_back(0);
  constant_eci.value.push_back(0);

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

  // jsonParser json;
  // json.put_array();
  // for (clexmonte::EventImpactInfo const &impact : prim_impact_info_list) {
  //   jsonParser tjson;
  //   to_json(impact.phenomenal_sites, tjson["phenomenal_sites"]);
  //   to_json(impact.required_update_neighborhood,
  //           tjson["required_update_neighborhood"]);
  //   json.push_back(tjson);
  // }
  // std::cout << json << std::endl;
}
