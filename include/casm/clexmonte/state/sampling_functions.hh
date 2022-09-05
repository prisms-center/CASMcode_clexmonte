#ifndef CASM_clexmonte_state_sampling_functions
#define CASM_clexmonte_state_sampling_functions

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

/// \brief Make temperature sampling function ("temperature")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_temperature_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make mol composition sampling function ("mol_composition")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_mol_composition_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make parametric composition sampling function ("param_composition")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_param_composition_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_corr_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make formation energy sampling function ("formation_energy")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_f(
    std::shared_ptr<SystemType> const &system);

/// \brief Make potential energy sampling function ("potential_energy")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_potential_energy_f(
    std::shared_ptr<SystemType> const &system);

// --- Inline definitions ---

/// \brief Make temperature sampling function ("temperature")
///
/// Requires:
/// - "temperature" is a scalar state condition
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_temperature_f(
    std::shared_ptr<SystemType> const &system) {
  return monte::StateSamplingFunction<Configuration>(
      "temperature", "Temperature (K)",
      1,  // number of components in "temperature",
      [](monte::State<Configuration> const &state) {
        return monte::reshaped(
            state.conditions.scalar_values.at("temperature"));
      });
}

/// \brief Make mol composition sampling function ("mol_composition")
///
/// Requires :
/// - `composition::CompositionConverter const &
///   get_composition_converter(SystemType &)`
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_mol_composition_f(
    std::shared_ptr<SystemType> const &system) {
  return monte::StateSamplingFunction<Configuration>(
      "mol_composition",
      "Number of each component (normalized per primitive cell)",
      get_composition_converter(*system).components(),  // component names
      [system](monte::State<Configuration> const &state) {
        Eigen::VectorXi const &occupation = get_occupation(state);
        return get_composition_calculator(*system).mean_num_each_component(
            occupation);
      });
}

/// \brief Make parametric composition sampling function ("param_composition")
///
/// Requires :
/// - `composition::CompositionConverter const &
///   get_composition_converter(SystemType &)`
/// - `composition::CompositionCalculator
///   get_composition_calculator(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_param_composition_f(
    std::shared_ptr<SystemType> const &system) {
  // name comp_x components "a", "b", ... for each independent composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(*system);
  std::vector<std::string> comp_x_components;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    comp_x_components.push_back(composition_converter.comp_var(i));
  }

  return monte::StateSamplingFunction<Configuration>(
      "param_composition", "Parametric composition",
      comp_x_components,  // component names
      [system](monte::State<Configuration> const &state) {
        composition::CompositionCalculator const &composition_calculator =
            get_composition_calculator(*system);
        composition::CompositionConverter const &composition_converter =
            get_composition_converter(*system);

        Eigen::VectorXi const &occupation = get_occupation(state);
        Eigen::VectorXd mol_composition =
            composition_calculator.mean_num_each_component(occupation);
        return composition_converter.param_composition(mol_composition);
      });
}

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
///
/// Requires:
/// - `ClexData &get_basis_set(SystemType &, std::string const &key)`
/// - `clexulator::ClusterExpansion &get_clex(SystemType &,
///   StateType const &, std::string const &key)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_corr_f(
    std::shared_ptr<SystemType> const &system) {
  // correlations size
  clexulator::Clexulator const &clexulator =
      *get_basis_set(*system, "formation_energy");
  Index corr_size = clexulator.corr_size();

  return monte::StateSamplingFunction<Configuration>(
      "formation_energy_corr",
      "Formation energy basis set correlations (normalized per primitive cell)",
      corr_size,  // number of components in "corr"
      [system](monte::State<Configuration> const &state) {
        clexulator::Correlations &correlations =
            get_clex(*system, state, "formation_energy")->correlations();
        auto const &extensive_corr = correlations.extensive();
        return correlations.intensive(extensive_corr);
      });
}

/// \brief Make formation energy sampling function ("formation_energy")
///
/// Requires:
/// - `clexulator::ClusterExpansion &get_clex(SystemType &,
///    StateType const &, std::string const &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_f(
    std::shared_ptr<SystemType> const &system) {
  return monte::StateSamplingFunction<Configuration>(
      "formation_energy",
      "Formation energy of the configuration (normalized per primitive cell)",
      1,  // number of components in "formation_energy"
      [system](monte::State<Configuration> const &state) {
        Eigen::VectorXd value(1);
        value(0) =
            get_clex(*system, state, "formation_energy")->intensive_value();
        return value;
      });
}

/// \brief Make potential energy sampling function ("potential_energy")
///
/// Requires:
/// - "potential_energy" is a scalar state property
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_potential_energy_f(
    std::shared_ptr<SystemType> const &system) {
  return monte::StateSamplingFunction<Configuration>(
      "potential_energy",
      "Potential energy of the state (normalized per primitive cell)",
      1,  // number of components in "potential_energy"
      [system](monte::State<Configuration> const &state) {
        return monte::reshaped(
            state.properties.scalar_values.at("potential_energy"));
      });
}

}  // namespace clexmonte
}  // namespace CASM

#endif
