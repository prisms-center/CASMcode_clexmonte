// #include "casm/monte/state/FixedConfigGenerator.hh"
#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"
#include "casm/monte/state/State.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class state_FixedConfigGeneratorTest : public test::ZrOTestSystem {};

TEST_F(state_FixedConfigGeneratorTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  EXPECT_EQ(get_basis_size(*system), 4);

  monte::ValueMap conditions =
      canonical::make_conditions(300.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});
  std::vector<State<Configuration>> finished_states;

  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system, T);
  for (Index i = 0; i < volume; ++i) {
    init_config.dof_values.occupation(2 * volume + i) = 1;
  }

  FixedConfigGenerator<Configuration> config_generator(init_config);
  for (Index j = 0; j < 10; ++j) {
    Configuration config = config_generator(conditions, finished_states);
    EXPECT_EQ(config.dof_values.occupation, init_config.dof_values.occupation);
  }
}
