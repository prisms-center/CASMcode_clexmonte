#ifndef CASM_clexmonte_monte_calculator_analysis_functions
#define CASM_clexmonte_monte_calculator_analysis_functions

#include "casm/clexmonte/definitions.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

namespace monte_calculator {

/// \brief Make heat capacity analysis function ("heat_capacity")
results_analysis_function_type make_heat_capacity_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make mol_composition susceptibility analysis function
/// ("mol_susc(A,B)")
results_analysis_function_type make_mol_susc_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make param_composition susceptibility analysis function
/// ("param_susc(a,b)")
results_analysis_function_type make_param_susc_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make mol_composition thermo-chemical susceptibility
///     analysis function ("mol_thermochem_susc(S,A)")
results_analysis_function_type make_mol_thermochem_susc_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make param_composition thermo-chemical susceptibility
///     analysis function ("param_thermochem_susc(S,a)")
results_analysis_function_type make_param_thermochem_susc_f(
    std::shared_ptr<MonteCalculator> const &calculation);

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
