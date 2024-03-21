#include "casm/clexmonte/semigrand_canonical/calculator_impl.hh"
#include "casm/clexmonte/state/make_conditions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/ValueMap.hh"
#include "casm/monte/run_management/ResultsAnalysisFunction.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

/// \brief Helper for making a conditions ValueMap for semi-grand
///     canonical Monte Carlo calculations
///
/// \param temperature The temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input.
/// \param param_chem_pot A map of axes names (for parametric composition)
///     to parametric chemical potential value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "param_chem_pot".
///
/// Example: Specifying "param_chem_pot"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                     // temperature (K)
///    composition_converter,     // composition converter
///    {{"a", -0.3}, {"b", 0.2}); // param_chem_pot values
/// \endcode
///
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> param_chem_pot) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["param_chem_pot"] =
      make_param_chem_pot(composition_converter, param_chem_pot);
  return conditions;
}

/// \brief Helper for making a conditions increment ValueMap for
///     semi-grand canonical Monte Carlo calculations
///
/// \param temperature The change in temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input.
/// \param param_chem_pot A map of axes names (for parametric composition)
///     to parametric chemical potential increment value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "param_chem_pot" (increment).
///
/// Example: Specifying "param_chem_pot"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                     // temperature (K)
///    composition_converter,     // composition converter
///    {{"a", 0.01}, {"b", 0.0}); // param_chem_pot increment values
/// \endcode
///
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> param_chem_pot) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["param_chem_pot"] =
      make_param_chem_pot_increment(composition_converter, param_chem_pot);
  return conditions;
}

template struct SemiGrandCanonical<std::mt19937_64>;

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM
