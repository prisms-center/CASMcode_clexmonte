#ifndef CASM_clexmonte_state_kinetic_sampling_functions
#define CASM_clexmonte_state_kinetic_sampling_functions

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/Correlations.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/state/StateSampler.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {

// ---
// These methods are used to construct sampling functions. They are templated
// so that they can be reused. The definition documentation should
// state interface requirements for the methods to be applicable and usable in
// a particular context.
//
// Example requirements are:
// - that a conditions `monte::ValueMap` contains scalar "temperature"
// - that the method `ClexData &get_clex(SystemType &,
//   StateType const &, std::string const &key)`
//   exists for template type `SystemType` (i.e. when
//   SystemType=clexmonte::System).
// ---

/// \brief Make KMC potential energy sampling function ("potential_energy")
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_kmc_potential_energy_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make center of mass squared displacement sampling function
/// ("R_squared_center")
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_R_squared_center_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make tracer squared displacement sampling function
/// ("R_squared_tracer")
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_R_squared_tracer_f(
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Make KMC potential energy sampling function ("potential_energy")
///
/// This is a canonical potential energy, equal to the formation energy,
/// but entirely calculated at the sampling time, instead of taken from
/// state.properties.
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_kmc_potential_energy_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return monte::StateSamplingFunction<Configuration>(
      "potential_energy",
      "Potential energy of the state (normalized per primitive cell)",
      {},  // scalar
      [calculation](monte::State<Configuration> const &state) {
        canonical::CanonicalPotential potential(calculation->system);
        potential.set(&state, calculation->conditions);
        double n_unitcells =
            get_transformation_matrix_to_supercell(state).determinant();
        return monte::reshaped(potential.extensive_value() / n_unitcells);
      });
}

/// \brief Make center of mass squared displacement sampling function
/// ("R_squared_center")
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_R_squared_center_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &components =
      get_composition_converter(*calculation->system).components();
  std::vector<std::string> dirs({"x", "y", "z"});
  std::vector<std::string> component_names;
  for (Index i = 0; i < components.size(); ++i) {
    for (Index j = i; j < components.size(); ++j) {
      for (Index alpha = 0; alpha < dirs.size(); ++alpha) {
        for (Index beta = alpha; beta < dirs.size(); ++beta) {
          component_names.push_back(components[i] + "," + components[j] + "," +
                                    dirs[alpha] + "," + dirs[beta]);
        }
      }
    }
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());
  return monte::StateSamplingFunction<Configuration>(
      "R_squared_center",
      R"(Samples \left(\sum_\zeta \Delta R^\zeta_i \right) \left(\sum_\zeta \Delta R^\zeta_j \right))",
      component_names,  // component names
      shape, [calculation](monte::State<Configuration> const &state) {
        auto const &components =
            get_composition_converter(*calculation->system).components();
        std::vector<std::string> dirs({"x", "y", "z"});
        std::vector<std::string> component_names;
        for (Index i = 0; i < components.size(); ++i) {
          for (Index j = i; j < components.size(); ++j) {
            for (Index alpha = 0; alpha < dirs.size(); ++alpha) {
              for (Index beta = alpha; beta < dirs.size(); ++beta) {
                component_names.push_back(components[i] + "," + components[j] +
                                          "," + dirs[alpha] + "," + dirs[beta]);
              }
            }
          }
        }
        Eigen::VectorXd v = Eigen::VectorXd::Zero(component_names.size());
        // TODO: ...
        return v;
      });
}

/// \brief Make tracer squared displacement sampling function
/// ("R_squared_tracer")
template <typename CalculationType>
monte::StateSamplingFunction<Configuration> make_R_squared_tracer_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &components =
      get_composition_converter(*calculation->system).components();
  std::vector<std::string> dirs({"x", "y", "z"});
  std::vector<std::string> component_names;
  for (Index i = 0; i < components.size(); ++i) {
    for (Index alpha = 0; alpha < dirs.size(); ++alpha) {
      for (Index beta = alpha; beta < dirs.size(); ++beta) {
        component_names.push_back(components[i] + "," + dirs[alpha] + "," +
                                  dirs[beta]);
      }
    }
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());
  return monte::StateSamplingFunction<Configuration>(
      "R_squared_tracer",
      R"(Samples \sum_\zeta \left(\Delta R^\zeta_{i,\alpha} \Delta R^\zeta_{i,\beta}\right)^2)",
      component_names,  // component names
      shape, [calculation](monte::State<Configuration> const &state) {
        auto const &components =
            get_composition_converter(*calculation->system).components();
        std::vector<std::string> dirs({"x", "y", "z"});
        std::vector<std::string> component_names;
        for (Index i = 0; i < components.size(); ++i) {
          for (Index alpha = 0; alpha < dirs.size(); ++alpha) {
            for (Index beta = alpha; beta < dirs.size(); ++beta) {
              component_names.push_back(components[i] + "," + dirs[alpha] +
                                        "," + dirs[beta]);
            }
          }
        }
        Eigen::VectorXd v = Eigen::VectorXd::Zero(component_names.size());
        // TODO: ...
        return v;
      });
}

}  // namespace clexmonte
}  // namespace CASM

#endif
