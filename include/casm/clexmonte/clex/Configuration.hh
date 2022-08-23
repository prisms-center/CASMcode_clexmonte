#ifndef CASM_clexmonte_clex_Configuration
#define CASM_clexmonte_clex_Configuration

#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"
#include "casm/monte/state/State.hh"

namespace CASM {
namespace clexmonte {

struct Configuration {
  Configuration(Eigen::Matrix3l const &_transformation_matrix_to_super,
                clexulator::ConfigDoFValues const &_dof_values)
      : transformation_matrix_to_super(_transformation_matrix_to_super),
        dof_values(_dof_values) {}

  Eigen::Matrix3l transformation_matrix_to_super;
  clexulator::ConfigDoFValues dof_values;
};

// --- The following are used to interface with CASM::monte methods ---

inline Eigen::Matrix3l const &get_transformation_matrix_to_super(
    Configuration const &configuration) {
  return configuration.transformation_matrix_to_super;
}

inline Eigen::VectorXi &get_occupation(Configuration &configuration) {
  return configuration.dof_values.occupation;
}

inline Eigen::VectorXi const &get_occupation(
    Configuration const &configuration) {
  return configuration.dof_values.occupation;
}

inline clexulator::ConfigDoFValues const &get_dof_values(
    Configuration const &configuration) {
  return configuration.dof_values;
}

inline clexulator::ConfigDoFValues &get_dof_values(
    Configuration &configuration) {
  return configuration.dof_values;
}

/// \brief Set calculator so it evaluates using `configuration`
inline void set(clexulator::ClusterExpansion &calculator,
                Configuration const &configuration) {
  calculator.set(&configuration.dof_values);
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::ClusterExpansion &calculator,
                monte::State<Configuration> const &state) {
  set(calculator, state.configuration);
}

/// \brief Set calculator so it evaluates using `configuration`
inline void set(clexulator::LocalClusterExpansion &calculator,
                Configuration const &configuration) {
  calculator.set(&configuration.dof_values);
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::LocalClusterExpansion &calculator,
                monte::State<Configuration> const &state) {
  set(calculator, state.configuration);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
