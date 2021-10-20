#ifndef CASM_clexmonte_system_sampling_functions
#define CASM_clexmonte_system_sampling_functions

#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/Correlations.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/state/StateSampler.hh"

namespace CASM {
namespace clexmonte {

/// \brief Make temperature sampling function ("temperature")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_temperature_f(
    std::shared_ptr<SystemType> const &system_data);

/// \brief Make mol composition sampling function ("comp_n")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_comp_n_f(
    std::shared_ptr<SystemType> const &system_data);

/// \brief Make parametric composition sampling function ("comp_x")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_comp_x_f(
    std::shared_ptr<SystemType> const &system_data);

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_corr_f(
    std::shared_ptr<SystemType> const &system_data);

/// \brief Make formation energy sampling function ("formation_energy")
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_f(
    std::shared_ptr<SystemType> const &system_data);

// --- Inline definitions ---

/// \brief Make temperature sampling function ("temperature")
///
/// Requires:
/// - "temperature" is a state condition
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_temperature_f(
    std::shared_ptr<SystemType> const &system_data) {
  return monte::StateSamplingFunction<Configuration>(
      "temperature", "Temperature (K)",
      1,  // number of components in "temperature",
      [=](monte::State<Configuration> const &state) {
        return state.conditions.at("temperature");
      });
}

/// \brief Make mol composition sampling function ("comp_n")
///
/// Requires :
/// - `composition::CompositionConverter const &
///   get_composition_converter(SystemType &)`
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_comp_n_f(
    std::shared_ptr<SystemType> const &system_data) {
  return monte::StateSamplingFunction<Configuration>(
      "comp_n", "Number of each component (normalized per primitive cell)",
      get_composition_converter(*system_data).components(),  // component names
      [&](monte::State<Configuration> const &state) {
        Eigen::VectorXi const &occupation = get_occupation(state.configuration);
        return get_composition_calculator(*system_data)
            .mean_num_each_component(occupation);
      });
}

/// \brief Make parametric composition sampling function ("comp_x")
///
/// Requires :
/// - `composition::CompositionConverter const &
///   get_composition_converter(SystemType &)`
/// - `composition::CompositionCalculator
///   get_composition_calculator(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_comp_x_f(
    std::shared_ptr<SystemType> const &system_data) {
  // name comp_x components "a", "b", ... for each independent composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(*system_data);
  std::vector<std::string> comp_x_components;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    comp_x_components.push_back(composition_converter.comp_var(i));
  }

  return monte::StateSamplingFunction<Configuration>(
      "comp_x", "Parametric composition",
      comp_x_components,  // component names
      [&](monte::State<Configuration> const &state) {
        composition::CompositionCalculator const &composition_calculator =
            get_composition_calculator(*system_data);
        composition::CompositionConverter const &composition_converter =
            get_composition_converter(*system_data);

        Eigen::VectorXi const &occupation = get_occupation(state.configuration);
        Eigen::VectorXd mol_composition =
            composition_calculator.mean_num_each_component(occupation);
        return composition_converter.param_composition(mol_composition);
      });
}

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
///
/// Requires:
/// - `ClexData &get_formation_energy_clex_data(SystemType &)`
/// - `Clex &get_formation_energy_clex(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_corr_f(
    std::shared_ptr<SystemType> const &system_data) {
  // correlations size
  clexulator::Clexulator const &clexulator =
      *get_formation_energy_clex_data(*system_data).clexulator;
  Index corr_size = clexulator.corr_size();

  return monte::StateSamplingFunction<Configuration>(
      "formation_energy_corr",
      "Formation energy basis set correlations (normalized per primitive cell)",
      corr_size,  // number of components in "corr"
      [&](monte::State<Configuration> const &state) {
        clexulator::Correlations &correlations =
            get_formation_energy_clex(*system_data, state).correlations();
        auto const &extensive_corr = correlations.extensive();
        return correlations.intensive(extensive_corr);
      });
}

/// \brief Make formation energy sampling function ("formation_energy")
///
/// Requires:
/// - `ClexData &get_formation_energy_clex_data(SystemType &)`
/// - `Clex &get_formation_energy_clex(SystemType &)`
template <typename SystemType>
monte::StateSamplingFunction<Configuration> make_formation_energy_f(
    std::shared_ptr<SystemType> const &system_data) {
  return monte::StateSamplingFunction<Configuration>(
      "formation_energy",
      "Formation energy of the configuration (normalized per primitive cell)",
      1,  // number of components in "formation_energy"
      [&](monte::State<Configuration> const &state) {
        Eigen::VectorXd value(1);
        value(1) =
            get_formation_energy_clex(*system_data, state).intensive_value();
        return value;
      });
}

}  // namespace clexmonte
}  // namespace CASM

#endif
