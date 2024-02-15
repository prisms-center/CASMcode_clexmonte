#include "ZrOTestSystem.hh"
#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/ValueMap.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class canonical_ConditionsTest : public test::ZrOTestSystem {};

TEST_F(canonical_ConditionsTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  monte::ValueMap conditions =
      canonical::make_conditions(300.0, system->composition_converter,
                                 {{"Zr", 2.0}, {"O", 1.0}, {"Va", 1.0}});

  EXPECT_EQ(conditions.scalar_values.size(), 1);
  EXPECT_EQ(conditions.scalar_values.at("temperature"), 300.0);
  EXPECT_EQ(conditions.vector_values.size(), 1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition").size(), 3);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(0), 2.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(1), 1.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(2), 1.0);
}

TEST_F(canonical_ConditionsTest, Test2) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  monte::ValueMap conditions = canonical::make_conditions_increment(
      10.0, system->composition_converter,
      {{"Zr", 0.0}, {"O", 0.1}, {"Va", -0.1}});

  EXPECT_EQ(conditions.scalar_values.size(), 1);
  EXPECT_EQ(conditions.scalar_values.at("temperature"), 10.0);
  EXPECT_EQ(conditions.vector_values.size(), 1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition").size(), 3);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(0), 0.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(1), -0.1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(2), 0.1);
}

TEST_F(canonical_ConditionsTest, Test3) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  monte::ValueMap conditions = canonical::make_conditions(
      300.0, system->composition_converter, {{"a", 0.5}});

  EXPECT_EQ(conditions.scalar_values.size(), 1);
  EXPECT_EQ(conditions.scalar_values.at("temperature"), 300.0);
  EXPECT_EQ(conditions.vector_values.size(), 1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition").size(), 3);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(0), 2.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(1), 1.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(2), 1.0);
}

TEST_F(canonical_ConditionsTest, Test4) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;

  monte::ValueMap conditions = canonical::make_conditions_increment(
      10.0, system->composition_converter, {{"a", 0.05}});

  EXPECT_EQ(conditions.scalar_values.size(), 1);
  EXPECT_EQ(conditions.scalar_values.at("temperature"), 10.0);
  EXPECT_EQ(conditions.vector_values.size(), 1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition").size(), 3);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(0), 0.0);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(1), -0.1);
  EXPECT_EQ(conditions.vector_values.at("mol_composition")(2), 0.1);
}
