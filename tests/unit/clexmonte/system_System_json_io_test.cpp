#include "KMCTestSystem.hh"
#include "gtest/gtest.h"

using namespace CASM;

/// NOTE:
/// - This test is designed to copy data to the same directory each time, so
///   that the Clexulators do not need to be re-compiled.
/// - To clear existing data, remove the directory:
//    CASM_test_projects/FCCBinaryVacancySystemJsonIOTest directory
class system_FCCBinaryVacancySystemJsonIOTest : public test::KMCTestSystem {
 protected:
  system_FCCBinaryVacancySystemJsonIOTest()
      : KMCTestSystem(
            "FCC_binary_vacancy", "FCCBinaryVacancySystemJsonIOTest",
            test::data_dir("clexmonte") / "kmc" / "system_template.json") {}
};

TEST_F(system_FCCBinaryVacancySystemJsonIOTest, Test1) {
  set_clex("formation_energy", "default", "formation_energy_eci.json");

  {
    fs::path event_relpath = fs::path("kmc_events") / "event.A_Va_1NN";
    set_local_basis_set("A_Va_1NN");
    set_event("A_Va_1NN", event_relpath / "kra_eci.json",
              event_relpath / "freq_eci.json");
  }

  {
    fs::path event_relpath = fs::path("kmc_events") / "event.B_Va_1NN";
    set_local_basis_set("B_Va_1NN");
    set_event("B_Va_1NN", event_relpath / "kra_eci.json",
              event_relpath / "freq_eci.json");
  }
  write_input();
  make_system();

  EXPECT_TRUE(system != nullptr);
  EXPECT_EQ(system->basis_sets.size(), 1);
  EXPECT_EQ(system->local_basis_sets.size(), 2);
  EXPECT_EQ(system->clex_data.size(), 1);
  EXPECT_EQ(system->local_multiclex_data.size(), 2);
  EXPECT_EQ(system->event_type_data.size(), 2);
}
