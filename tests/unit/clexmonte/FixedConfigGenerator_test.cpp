// #include "casm/monte/state/FixedConfigGenerator.hh"
#include "ZrOTestSystem.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class FixedConfigGeneratorTest : public test::ZrOTestSystem {};

TEST_F(FixedConfigGeneratorTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  EXPECT_EQ(system_data->shared_prim->basis().size(), 4);

  // VectorValueMap conditions = canonical::make_conditions();
  //
  // FixedConfigGenerator<clexmonte::Configuration> generator(config);
}
