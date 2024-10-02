#ifndef CASM_clexmonte_state_Configuration
#define CASM_clexmonte_state_Configuration

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/LocalOrbitCompositionCalculator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

typedef config::Configuration Configuration;

// struct Configuration {
//   Configuration(Eigen::Matrix3l const &_transformation_matrix_to_super,
//                 clexulator::ConfigDoFValues const &_dof_values)
//       : transformation_matrix_to_super(_transformation_matrix_to_super),
//         dof_values(_dof_values) {}
//
//   Eigen::Matrix3l transformation_matrix_to_super;
//   clexulator::ConfigDoFValues dof_values;
// };

// --- The following are used to interface with CASM::monte methods ---

inline Eigen::Matrix3l const &get_transformation_matrix_to_super(
    state_type const &state) {
  return state.configuration.supercell->superlattice
      .transformation_matrix_to_super();
}

inline Eigen::VectorXi &get_occupation(state_type &state) {
  return state.configuration.dof_values.occupation;
}

inline Eigen::VectorXi const &get_occupation(state_type const &state) {
  return state.configuration.dof_values.occupation;
}

inline clexulator::ConfigDoFValues const &get_dof_values(
    state_type const &state) {
  return state.configuration.dof_values;
}

inline clexulator::ConfigDoFValues &get_dof_values(state_type &state) {
  return state.configuration.dof_values;
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::ClusterExpansion &calculator,
                state_type const &state) {
  calculator.set(&get_dof_values(state));
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::MultiClusterExpansion &calculator,
                state_type const &state) {
  calculator.set(&get_dof_values(state));
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::LocalClusterExpansion &calculator,
                state_type const &state) {
  calculator.set(&get_dof_values(state));
}

/// \brief Set calculator so it evaluates using `state`
inline void set(clexulator::MultiLocalClusterExpansion &calculator,
                state_type const &state) {
  calculator.set(&get_dof_values(state));
}

/// \brief Set calculator so it evaluates using `state`
inline void set(LocalOrbitCompositionCalculator &calculator,
                state_type const &state) {
  calculator.set(&get_dof_values(state));
}

}  // namespace clexmonte

namespace monte {
using clexmonte::get_occupation;
using clexmonte::get_transformation_matrix_to_super;
}  // namespace monte
}  // namespace CASM

#endif
