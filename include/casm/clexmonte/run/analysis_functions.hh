#ifndef CASM_clexmonte_run_analysis_functions
#define CASM_clexmonte_run_analysis_functions

#include "casm/clexmonte/run/covariance_functions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {

// ---
// These methods are used to construct results analysis functions. They are
// templated so that they can be reused. The definition documentation should
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

/// \brief Make heat capacity analysis function ("heat_capacity")
template <typename CalculationType>
results_analysis_function_type make_heat_capacity_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make mol_composition susceptibility analysis function
/// ("mol_susc(A,B)")
template <typename CalculationType>
results_analysis_function_type make_mol_susc_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make param_composition susceptibility analysis function
/// ("param_susc(a,b)")
template <typename CalculationType>
results_analysis_function_type make_param_susc_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make mol_composition thermo-chemical susceptibility
///     analysis function ("mol_thermochem_susc(S,A)")
template <typename CalculationType>
results_analysis_function_type make_mol_thermochem_susc_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make param_composition thermo-chemical susceptibility
///     analysis function ("param_thermochem_susc(S,a)")
template <typename CalculationType>
results_analysis_function_type make_param_thermochem_susc_f(
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Calculates `(kB * temperature * temperature) / n_unitcells`
template <typename CalculationType>
std::function<double()> make_heat_capacity_normalization_constant_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return [=]() -> double {
    // validate temperature
    auto const &state = *calculation->state;
    auto const &conditions = state.conditions;
    Index n_unitcells = get_transformation_matrix_to_super(state).determinant();
    if (!conditions.scalar_values.count("temperature")) {
      std::stringstream msg;
      msg << "Results analysis error: heat_capacity requires temperature "
             "condition";
      throw std::runtime_error(msg.str());
    }
    double temperature = conditions.scalar_values.at("temperature");

    // calculate
    return (CASM::KB * temperature * temperature) / n_unitcells;
  };
}

/// \brief Make heat capacity analysis function ("heat_capacity")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename CalculationType>
results_analysis_function_type make_heat_capacity_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return make_variance_f(
      "heat_capacity",
      "Heat capacity (per unit cell) = "
      "var(potential_energy_per_unitcell)*n_unitcells/(kB*T*T)",
      "potential_energy", {"0"}, {},
      make_heat_capacity_normalization_constant_f<CalculationType>(
          calculation));
}

/// \brief Calculates `(kB * temperature) / n_unitcells`
template <typename CalculationType>
std::function<double()> make_susc_normalization_constant_f(
    std::shared_ptr<CalculationType> const &calculation, std::string name) {
  return [=]() -> double {
    // validate temperature
    auto const &state = *calculation->state;
    auto const &conditions = state.conditions;
    Index n_unitcells = get_transformation_matrix_to_super(state).determinant();
    if (!conditions.scalar_values.count("temperature")) {
      std::stringstream msg;
      msg << "Results analysis error: " << name
          << " requires temperature condition";
      throw std::runtime_error(msg.str());
    }
    double temperature = conditions.scalar_values.at("temperature");

    // calculate
    return (CASM::KB * temperature) / n_unitcells;
  };
}

/// \brief Make mol_composition susceptibility analysis function
/// ("mol_susc(A,B)")
///
/// Notes:
/// - Requires sampling "mol_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename CalculationType>
results_analysis_function_type make_mol_susc_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  auto const &component_names = get_composition_converter(system).components();
  return make_covariance_f(
      "mol_susc",
      "Chemical susceptibility (per unit cell) = "
      "cov(mol_composition_i, mol_composition_j)*n_unitcells/(kB*T)",
      "mol_composition", "mol_composition", component_names, component_names,
      make_susc_normalization_constant_f(calculation, "mol_susc"));
}

/// \brief Make param_composition susceptibility analysis function
/// ("param_susc(a,b)")
///
/// Notes:
/// - Requires sampling "param_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename CalculationType>
results_analysis_function_type make_param_susc_f(
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
  return make_covariance_f(
      "param_susc",
      "Chemical susceptibility (per unit cell) = "
      "cov(param_composition_i, param_composition_j)*n_unitcells/(kB*T)",
      "param_composition", "param_composition", component_names,
      component_names,
      make_susc_normalization_constant_f(calculation, "param_susc"));
}

/// \brief Make mol_composition thermo-chemical susceptibility
///     analysis function ("mol_thermochem_susc(S,A)")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires sampling "mol_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename CalculationType>
results_analysis_function_type make_mol_thermochem_susc_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  std::vector<std::string> first_component_names = {"S"};

  auto const &second_component_names =
      get_composition_converter(system).components();

  return make_covariance_f(
      "mol_thermochem_susc",
      "Thermo-chemical susceptibility (per unit cell) = "
      "cov(potential_energy, mol_composition)*n_unitcells/(kB*T)",
      "potential_energy", "mol_composition", first_component_names,
      second_component_names,
      make_susc_normalization_constant_f(calculation, "mol_thermochem_susc"));
}

/// \brief Make param_composition thermo-chemical susceptibility
///     analysis function ("param_thermochem_susc(S,a)")
///
/// Notes:
/// - Requires sampling "potential_energy" (as per unit cell energy)
/// - Requires sampling "param_composition"
/// - Requires scalar condition "temperature"
/// - Requires result "initial_state"
template <typename CalculationType>
results_analysis_function_type make_param_thermochem_susc_f(
    std::shared_ptr<CalculationType> const &calculation) {
  auto const &system = *calculation->system;
  std::vector<std::string> first_component_names = {"S"};

  // name param_composition components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(system);
  std::vector<std::string> second_component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    second_component_names.push_back(composition_converter.comp_var(i));
  }
  return make_covariance_f(
      "param_thermochem_susc",
      "Thermo-chemical susceptibility (per unit cell) = "
      "cov(potential_energy, param_composition)*n_unitcells/(kB*T)",
      "potential_energy", "param_composition", first_component_names,
      second_component_names,
      make_susc_normalization_constant_f(calculation, "param_thermochem_susc"));
}

}  // namespace clexmonte
}  // namespace CASM

#endif
