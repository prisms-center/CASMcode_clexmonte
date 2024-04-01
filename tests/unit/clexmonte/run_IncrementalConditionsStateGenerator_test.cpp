#include "ZrOTestSystem.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/run/FixedConfigGenerator.hh"
#include "casm/clexmonte/run/IncrementalConditionsStateGenerator.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/monte/run_management/State.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

using namespace test;

class run_IncrementalConditionsStateGeneratorTest : public test::ZrOTestSystem {
};

TEST_F(run_IncrementalConditionsStateGeneratorTest, Test1) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;
  typedef FixedConfigGenerator fixed_config_generator_type;
  typedef IncrementalConditionsStateGenerator incremental_state_generator_type;

  EXPECT_EQ(get_basis_size(*system), 4);

  ValueMap init_conditions =
      canonical::make_conditions(300.0, get_composition_converter(*system),
                                 {{"Zr", 2.0}, {"O", 0.2}, {"Va", 1.8}});
  ValueMap conditions_increment = canonical::make_conditions_increment(
      10.0, get_composition_converter(*system),
      {{"Zr", 0.0}, {"O", 0.0}, {"Va", 0.0}});
  Index n_states = 11;
  bool dependent_runs = false;

  // init_config (for config_generator)
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system, T);
  for (Index i = 0; i < volume; ++i) {
    init_config.dof_values.occupation(2 * volume + i) = 1;
  }
  // config_generator
  std::unique_ptr<config_generator_type> config_generator =
      notstd::make_unique<fixed_config_generator_type>(init_config);

  // modifiers
  std::vector<StateModifyingFunction> modifiers;

  RunDataOutputParams output_params;
  output_params.output_dir = test_dir / "output";
  incremental_state_generator_type state_generator(
      system, output_params, std::move(config_generator), init_conditions,
      conditions_increment, n_states, dependent_runs, modifiers);

  while (!state_generator.is_complete()) {
    state_type state = state_generator.next_state();
    EXPECT_EQ(get_occupation(state), init_config.dof_values.occupation);
    EXPECT_TRUE(
        CASM::almost_equal(state.conditions.scalar_values.at("temperature"),
                           300.0 + 10.0 * state_generator.n_completed_runs()));
    RunData run_data;
    run_data.initial_state = state;
    run_data.final_state = state;
    run_data.conditions = state.conditions;
    run_data.transformation_matrix_to_super = T;
    run_data.n_unitcells = T.determinant();
    state_generator.push_back(run_data);
  }

  state_generator.write_completed_runs();

  EXPECT_TRUE(fs::exists(output_params.output_dir / "completed_runs.json"));
  {
    // test `read_completed_runs()`
    incremental_state_generator_type state_generator_read_test(
        system, output_params,
        notstd::make_unique<fixed_config_generator_type>(init_config),
        init_conditions, conditions_increment, n_states, dependent_runs,
        modifiers);
    state_generator_read_test.read_completed_runs();
    EXPECT_EQ(state_generator_read_test.n_completed_runs(), n_states);
  }
  fs::remove(output_params.output_dir / "completed_runs.json");
  fs::remove(output_params.output_dir);
}

TEST_F(run_IncrementalConditionsStateGeneratorTest, Test2) {
  using namespace CASM;
  using namespace CASM::monte;
  using namespace CASM::clexmonte;
  typedef FixedConfigGenerator fixed_config_generator_type;
  typedef IncrementalConditionsStateGenerator incremental_state_generator_type;

  EXPECT_EQ(get_basis_size(*system), 4);

  ValueMap init_conditions =
      canonical::make_conditions(300.0, get_composition_converter(*system),
                                 {{"Zr", 2.0}, {"O", 0.2}, {"Va", 1.8}});
  ValueMap conditions_increment = canonical::make_conditions_increment(
      0.0, get_composition_converter(*system),
      {{"Zr", 0.0}, {"O", 0.2}, {"Va", -0.2}});
  Index n_states = 9;
  bool dependent_runs = false;

  // init_config (for config_generator)
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
  Index volume = T.determinant();
  Configuration init_config = make_default_configuration(*system, T);
  for (Index i = 0; i < volume; ++i) {
    init_config.dof_values.occupation(2 * volume + i) = 1;
  }
  // config_generator
  std::unique_ptr<config_generator_type> config_generator =
      notstd::make_unique<fixed_config_generator_type>(init_config);

  // modifiers
  std::vector<StateModifyingFunction> modifiers;

  RunDataOutputParams output_params;
  incremental_state_generator_type state_generator(
      system, output_params, std::move(config_generator), init_conditions,
      conditions_increment, n_states, dependent_runs, modifiers);

  while (!state_generator.is_complete()) {
    state_type state = state_generator.next_state();
    EXPECT_EQ(get_occupation(state), init_config.dof_values.occupation);
    EXPECT_TRUE(almost_equal(
        state.conditions.vector_values.at("mol_composition"),
        init_conditions.vector_values.at("mol_composition") +
            conditions_increment.vector_values.at("mol_composition") *
                state_generator.n_completed_runs()));
    RunData run_data;
    run_data.initial_state = state;
    run_data.final_state = state;
    run_data.conditions = state.conditions;
    run_data.transformation_matrix_to_super = T;
    run_data.n_unitcells = T.determinant();
    state_generator.push_back(run_data);
  }
}
