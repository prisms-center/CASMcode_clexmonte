#ifndef CASM_clexmonte_monte_calculator_sampling_functions
#define CASM_clexmonte_monte_calculator_sampling_functions

#include "casm/clexmonte/definitions.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

namespace monte_calculator {

/// \brief Make temperature sampling function ("temperature")
///
/// Requires:
/// - "temperature" is a scalar state condition
state_sampling_function_type make_temperature_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make mol composition sampling function ("mol_composition")
state_sampling_function_type make_mol_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make parametric composition sampling function ("param_composition")
state_sampling_function_type make_param_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make parametric chemical potential sampling function
///     ("param_chem_pot")
///
/// Requires:
/// - "param_chem_pot" is a vector state condition
state_sampling_function_type make_param_chem_pot_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make correlations sampling function ("corr.<key>")
state_sampling_function_type make_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make cluster expansion value sampling function ("clex.<key>")
state_sampling_function_type make_clex_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make multi-cluster expansion value sampling function
/// ("multiclex.<key>")
state_sampling_function_type make_multiclex_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make formation energy correlations sampling function
///     ("clex.<key>.sparse_corr")
state_sampling_function_type make_clex_sparse_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make formation energy correlations sampling function
///     ("multiclex.<key>.sparse_corr")
state_sampling_function_type make_multiclex_sparse_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make potential energy sampling function (allows user-specified
/// <label>)
state_sampling_function_type make_potential_energy_f(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string label = "generalized_enthalpy",
    std::string desc =
        "Genaralized enthalpy of the state (normalized per primitive cell)");

/// \brief Make order parameter sampling function ("order_parameter.<key>")
state_sampling_function_type make_order_parameter_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make order parameter magnitudes by subspace sampling function
/// ("order_parameter.<key>.subspace_magnitudes")
state_sampling_function_type make_subspace_order_parameter_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key);

/// \brief Make configuration sampling function
json_state_sampling_function_type make_config_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make sampling functions that all methods can use
std::vector<state_sampling_function_type> common_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculations,
    std::string potential_energy_name = "generalized_enthalpy",
    std::string potential_energy_desc =
        "Genaralized enthalpy of the state (normalized per primitive cell)");

/// \brief Make json sampling functions that all methods can use
std::vector<json_state_sampling_function_type> common_json_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculations);

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
