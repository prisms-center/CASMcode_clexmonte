#include "KMCTestSystem.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "gtest/gtest.h"

// write neighborhoods to json:
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancySystemJsonIOTest directory
class kmc_FCCBinaryVacancySystemJsonIOTest : public KMCTestSystem {
 protected:
  kmc_FCCBinaryVacancySystemJsonIOTest()
      : KMCTestSystem(
            "FCC_binary_vacancy", "FCCBinaryVacancySystemJsonIOTest",
            test::data_dir("clexmonte") / "kmc" / "system_template.json") {}
};

TEST_F(kmc_FCCBinaryVacancySystemJsonIOTest, Test1) {
  set_clex("formation_energy", "default", "formation_energy_eci.json");

  {
    fs::path event_relpath = fs::path("events") / "event.A_Va_1NN";
    set_local_basis_set("A_Va_1NN");
    set_event("A_Va_1NN", event_relpath / "kra_eci.json",
              event_relpath / "freq_eci.json");
  }

  {
    fs::path event_relpath = fs::path("events") / "event.B_Va_1NN";
    set_local_basis_set("B_Va_1NN");
    set_event("B_Va_1NN", event_relpath / "kra_eci.json",
              event_relpath / "freq_eci.json");
  }
  write_input();

  ParentInputParser parser(json["kwargs"]);
  auto subparser = parser.subparse<clexmonte::System>("system");
  std::runtime_error error_if_invalid{"Error reading KMC System JSON input"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  EXPECT_TRUE(subparser->valid());

  clexmonte::System &system = *subparser->value;
  EXPECT_EQ(system.basis_sets.size(), 1);
  EXPECT_EQ(system.local_basis_sets.size(), 2);
  EXPECT_EQ(system.clex_data.size(), 1);
  EXPECT_EQ(system.local_multiclex_data.size(), 2);
  EXPECT_EQ(system.event_type_data.size(), 2);

  // jsonParser njson;
  // njson["formation_energy"] = get_required_update_neighborhood(
  //     system, system.clex_data.at("formation_energy"));
  // for (auto const &pair : system.event_type_data) {
  //   auto const &local_multiclex_data =
  //       system.local_multiclex_data.at(pair.second.local_multiclex_name);
  //   auto n_equivalents = pair.second.events.size();
  //   for (Index equivalent_index = 0; equivalent_index < n_equivalents;
  //        ++equivalent_index) {
  //     njson["local_clex"][pair.first][std::to_string(equivalent_index)] =
  //         get_required_update_neighborhood(system, local_multiclex_data,
  //                                          equivalent_index);
  //   }
  // }
  // std::cout << "update neighborhoods:" << std::endl;
  // std::cout << njson << std::endl;
}
