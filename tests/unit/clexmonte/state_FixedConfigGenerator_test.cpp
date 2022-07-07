// #include "casm/monte/state/FixedConfigGenerator.hh"
#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/conditions.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"
#include "casm/monte/state/State.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class FixedConfigGeneratorTest : public test::ZrOTestSystem {};

TEST_F(FixedConfigGeneratorTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  EXPECT_EQ(system_data->shared_prim->basis().size(), 4);

  monte::VectorValueMap conditions =
      canonical::make_conditions(300.0, system_data->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});
  std::vector<State<Configuration>> finished_states;

  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system_data, T);
  for (Index i = 0; i < volume; ++i) {
    get_occupation(init_config)(2 * volume + i) = 1;
  }

  FixedConfigGenerator<Configuration> config_generator(init_config);
  for (Index j = 0; j < 10; ++j) {
    Configuration config = config_generator(conditions, finished_states);
    EXPECT_EQ(get_occupation(config), get_occupation(init_config));
  }
}
