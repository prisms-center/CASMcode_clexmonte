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
template <typename CalculationType>
state_sampling_function_type make_temperature_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make mol composition sampling function ("mol_composition")
template <typename CalculationType>
state_sampling_function_type make_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make parametric composition sampling function ("param_composition")
template <typename CalculationType>
state_sampling_function_type make_param_composition_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make parametric chemical potential sampling function
/// ("param_chem_pot")
template <typename CalculationType>
state_sampling_function_type make_param_chem_pot_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
template <typename CalculationType>
state_sampling_function_type make_formation_energy_corr_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make formation energy sampling function ("formation_energy")
template <typename CalculationType>
state_sampling_function_type make_formation_energy_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make potential energy sampling function ("potential_energy")
template <typename CalculationType>
state_sampling_function_type make_potential_energy_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make order parameter sampling function ("order_parameter_X")
template <typename CalculationType>
void make_order_parameter_f(
    std::vector<state_sampling_function_type> &functions,
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Make temperature sampling function ("temperature")
///
/// Requires:
/// - "temperature" is a scalar state condition
template <typename CalculationType>
state_sampling_function_type make_temperature_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return state_sampling_function_type(
      "temperature", "Temperature (K)", {},  // scalar,
      [calculation]() {
        return monte::reshaped(calculation->conditions->temperature);
      });
}

/// \brief Make mol composition sampling function ("mol_composition")
///
/// Requires :
/// - `composition::CompositionConverter const &
///   get_composition_converter(SystemType &)`
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename CalculationType>
state_sampling_function_type make_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  auto const &components = get_composition_converter(system).components();
  std::vector<Index> shape;
  shape.push_back(components.size());
  return state_sampling_function_type(
      "mol_composition",
      "Number of each component (normalized per primitive cell)",
      components,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system;
        auto const &state = *calculation->state;
        Eigen::VectorXi const &occupation = get_occupation(state);
        return get_composition_calculator(system).mean_num_each_component(
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
template <typename CalculationType>
state_sampling_function_type make_param_composition_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  // name param_composition components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(system);
  std::vector<std::string> component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    component_names.push_back(composition_converter.comp_var(i));
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "param_composition", "Parametric composition",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system;
        auto const &state = *calculation->state;
        composition::CompositionCalculator const &composition_calculator =
            get_composition_calculator(system);
        composition::CompositionConverter const &composition_converter =
            get_composition_converter(system);

        Eigen::VectorXi const &occupation = get_occupation(state);
        Eigen::VectorXd mol_composition =
            composition_calculator.mean_num_each_component(occupation);
        return composition_converter.param_composition(mol_composition);
      });
}

/// \brief Make parametric chemical potential sampling function
/// ("param_chem_pot")
template <typename CalculationType>
state_sampling_function_type make_param_chem_pot_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  // name param_chem_pot components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(system);
  std::vector<std::string> component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    component_names.push_back(composition_converter.comp_var(i));
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "param_chem_pot",
      "Chemical potential conjugate to parametric composition axes",
      component_names,  // component names
      shape, [calculation]() {
        return calculation->conditions->param_chem_pot.value();
      });
}

/// \brief Make formation energy correlations sampling function
///     ("formation_energy_corr")
///
/// Requires:
/// - `ClexData &get_basis_set(SystemType &, std::string const &key)`
/// - `clexulator::ClusterExpansion &get_clex(SystemType &,
///   StateType const &, std::string const &key)`
template <typename CalculationType>
state_sampling_function_type make_formation_energy_corr_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  // correlations size
  clexulator::Clexulator const &clexulator =
      *get_basis_set(system, "formation_energy");
  Index corr_size = clexulator.corr_size();
  std::vector<Index> shape;
  shape.push_back(corr_size);

  return state_sampling_function_type(
      "formation_energy_corr",
      "Formation energy basis set correlations (normalized per primitive cell)",
      shape, [calculation]() {
        auto &system = *calculation->system;
        auto const &state = *calculation->state;
        clexulator::Correlations &correlations =
            *get_corr(system, state, "formation_energy");
        auto const &extensive_corr = correlations.extensive();
        return correlations.intensive(extensive_corr);
      });
}

/// \brief Make formation energy sampling function ("formation_energy")
///
/// Requires:
/// - `clexulator::ClusterExpansion &get_clex(SystemType &,
///    StateType const &, std::string const &)`
template <typename CalculationType>
state_sampling_function_type make_formation_energy_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return state_sampling_function_type(
      "formation_energy",
      "Formation energy of the configuration (normalized per primitive cell)",
      {},  // scalar
      [calculation]() {
        auto &system = *calculation->system;
        auto const &state = *calculation->state;
        Eigen::VectorXd value(1);
        value(0) =
            get_clex(system, state, "formation_energy")->intensive_value();
        return value;
      });
}

/// \brief Make potential energy sampling function ("potential_energy")
///
/// Notes:
/// - This version reads from state.properties, so it works
///   for methods such as `Canonical` and `SemiGrandCanonical`
///   which keep the potential_energy updated, but not `Kinetic`.
///   For `Kinetic` use `make_canonical_potential_energy_f`.
///
/// Requires:
/// - "potential_energy" is a scalar state property
template <typename CalculationType>
state_sampling_function_type make_potential_energy_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return state_sampling_function_type(
      "potential_energy",
      "Potential energy of the state (normalized per primitive cell)",
      {},  // scalar
      [calculation]() {
        auto const &state = *calculation->state;
        return monte::reshaped(
            state.properties.scalar_values.at("potential_energy"));
      });
}

/// \brief Make order parameter sampling function ("order_parameter_<key>")
///
/// Creates one "order_parameter_<key>" function for each DoFSpace in
/// the `calculation->system->dof_spaces` map.
template <typename CalculationType>
void make_order_parameter_f(
    std::vector<state_sampling_function_type> &functions,
    std::shared_ptr<CalculationType> const &calculation) {
  for (auto const &pair : calculation->system->dof_spaces) {
    std::string key = pair.first;
    clexulator::DoFSpace const &dof_space = pair.second;
    std::string name = "order_parameter_" + key;
    std::string desc = "Order parameters";

    functions.push_back(state_sampling_function_type(
        name, desc, {dof_space.subspace_dim},  // vector size
        [calculation, key]() {
          return get_order_parameter(*calculation->system, *calculation->state,
                                     key)
              ->value();
        }));
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
