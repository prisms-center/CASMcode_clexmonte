#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/conditions.hh"
#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/monte/state/FixedConfigGenerator.hh"
#include "casm/monte/state/IncrementalConditionsStateGenerator.hh"
#include "casm/monte/state/State.hh"
#include "casm/monte/state/StateSampler.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class IncrementalConditionsStateGeneratorTest : public test::ZrOTestSystem {};

TEST_F(IncrementalConditionsStateGeneratorTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;
  typedef Configuration config_type;
  typedef State<config_type> state_type;
  typedef ConfigGenerator<config_type, state_type> config_generator_type;
  typedef FixedConfigGenerator<config_type> fixed_config_generator_type;
  typedef IncrementalConditionsStateGenerator<config_type>
      incremental_state_generator_type;

  EXPECT_EQ(system_data->shared_prim->basis().size(), 4);

  VectorValueMap init_conditions =
      canonical::make_conditions(300.0, system_data->composition_converter,
                                 {{"Zr", 2.0}, {"O", 0.2}, {"Va", 1.8}});
  VectorValueMap conditions_increment = canonical::make_conditions_increment(
      10.0, system_data->composition_converter,
      {{"Zr", 0.0}, {"O", 0.0}, {"Va", 0.0}});
  Index n_states = 11;
  bool dependent_runs = false;

  // init_config (for config_generator)
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system_data, T);
  for (Index i = 0; i < volume; ++i) {
    init_config.dof_values.occupation(2 * volume + i) = 1;
  }
  // config_generator
  std::unique_ptr<config_generator_type> config_generator =
      notstd::make_unique<fixed_config_generator_type>(init_config);

  // dependent_conditions
  StateSamplingFunctionMap<config_type> dependent_conditions;

  std::vector<state_type> final_states;
  incremental_state_generator_type state_generator(
      std::move(config_generator), init_conditions, conditions_increment,
      n_states, dependent_runs, dependent_conditions);

  while (!state_generator.is_complete(final_states)) {
    state_type state = state_generator.next_state(final_states);
    Configuration const &config = state.configuration;
    EXPECT_EQ(config.dof_values.occupation, init_config.dof_values.occupation);
    EXPECT_TRUE(CASM::almost_equal(state.conditions.at("temperature")(0),
                                   300.0 + 10.0 * final_states.size()));
    final_states.push_back(state);
  }
}

TEST_F(IncrementalConditionsStateGeneratorTest, Test2) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;
  typedef Configuration config_type;
  typedef State<config_type> state_type;
  typedef ConfigGenerator<config_type, state_type> config_generator_type;
  typedef FixedConfigGenerator<config_type> fixed_config_generator_type;
  typedef IncrementalConditionsStateGenerator<config_type>
      incremental_state_generator_type;

  EXPECT_EQ(system_data->shared_prim->basis().size(), 4);

  VectorValueMap init_conditions =
      canonical::make_conditions(300.0, system_data->composition_converter,
                                 {{"Zr", 2.0}, {"O", 0.2}, {"Va", 1.8}});
  VectorValueMap conditions_increment = canonical::make_conditions_increment(
      0.0, system_data->composition_converter,
      {{"Zr", 0.0}, {"O", 0.2}, {"Va", -0.2}});
  Index n_states = 9;
  bool dependent_runs = false;

  // init_config (for config_generator)
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system_data, T);
  for (Index i = 0; i < volume; ++i) {
    init_config.dof_values.occupation(2 * volume + i) = 1;
  }
  // config_generator
  std::unique_ptr<config_generator_type> config_generator =
      notstd::make_unique<fixed_config_generator_type>(init_config);

  // dependent_conditions
  StateSamplingFunctionMap<config_type> dependent_conditions;

  std::vector<state_type> final_states;
  incremental_state_generator_type state_generator(
      std::move(config_generator), init_conditions, conditions_increment,
      n_states, dependent_runs, dependent_conditions);

  while (!state_generator.is_complete(final_states)) {
    state_type state = state_generator.next_state(final_states);
    Configuration const &config = state.configuration;
    EXPECT_EQ(config.dof_values.occupation, init_config.dof_values.occupation);
    EXPECT_TRUE(almost_equal(
        state.conditions.at("mol_composition"),
        init_conditions.at("mol_composition") +
            conditions_increment.at("mol_composition") * final_states.size()));
    final_states.push_back(state);
  }
}
