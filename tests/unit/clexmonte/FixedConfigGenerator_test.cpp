// #include "casm/monte/state/FixedConfigGenerator.hh"
#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/conditions.hh"
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

  monte::VectorValueMap conditions =
      canonical::make_conditions(300.0, system_data->composition_converter,
                                 {{"Zr", 0.0}, {"O", 1.0}, {"Va", 1.0}});

  jsonParser json;
  to_json(conditions, json["conditions"], jsonParser::as_flattest());
  std::cout << json << std::endl;

  // FixedConfigGenerator<clexmonte::Configuration> generator(config);
}
