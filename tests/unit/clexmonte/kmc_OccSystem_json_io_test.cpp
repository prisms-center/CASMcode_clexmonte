#include "KMCTestSystem.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/io/json/OccSystem_json_io.hh"
#include "gtest/gtest.h"

// write neighborhoods to json:
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancyOccSystemJsonIOTest directory
class FCCBinaryVacancyOccSystemJsonIOTest : public KMCTestSystem {
 protected:
  FCCBinaryVacancyOccSystemJsonIOTest()
      : KMCTestSystem("FCC_binary_vacancy",
                      "FCCBinaryVacancyOccSystemJsonIOTest",
                      test::data_dir("clexmonte") / "kmc" / "system_1.json") {}
};

TEST_F(FCCBinaryVacancyOccSystemJsonIOTest, Test1) {
  set_clex("formation_energy", "default", "formation_energy_eci.json");

  {
    fs::path event_relpath = fs::path("events") / "event.A_Va_1NN";
    set_local_clex("A_Va_1NN_kra", "A_Va_1NN", event_relpath / "kra_eci.json");
    set_local_clex("A_Va_1NN_freq", "A_Va_1NN",
                   event_relpath / "freq_eci.json");
    set_event("A_Va_1NN", "A_Va_1NN_kra", "A_Va_1NN_freq");
  }

  {
    fs::path event_relpath = fs::path("events") / "event.B_Va_1NN";
    set_local_clex("B_Va_1NN_kra", "B_Va_1NN", event_relpath / "kra_eci.json");
    set_local_clex("B_Va_1NN_freq", "B_Va_1NN",
                   event_relpath / "freq_eci.json");
    set_event("B_Va_1NN", "B_Va_1NN_kra", "B_Va_1NN_freq");
  }

  write_input();

  ParentInputParser parser(json["kwargs"]);
  auto subparser = parser.subparse<clexmonte::OccSystem>("system");
  std::runtime_error error_if_invalid{"Error reading KMC OccSystem JSON input"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  EXPECT_TRUE(subparser->valid());

  clexmonte::OccSystem &system = *subparser->value;
  jsonParser njson;
  njson["formation_energy"] =
      system.formation_energy_clex_data.required_update_neighborhood();
  for (auto const &pair : system.local_clex_data) {
    for (Index equivalent_index = 0;
         equivalent_index < pair.second.clexulator->size();
         ++equivalent_index) {
      njson["local_clex"][pair.first][std::to_string(equivalent_index)] =
          pair.second.required_update_neighborhood(equivalent_index);
    }
  }
  std::cout << "update neighborhoods:" << std::endl;
  std::cout << njson << std::endl;
}
