#include "casm/clexmonte/canonical/conditions.hh"

#include "casm/clexmonte/system/make_conditions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Helper for making a conditions VectorValueMap for canonical Monte
///     Carlo calculations
///
/// \param temperature The temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate `comp` input and convert to  `comp_n` (mol_composition).
/// \param comp A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
///
/// \returns VectorValueMap which contains "temperature" and "mol_composition".
///
/// Example: Specifying "mol_composition"
/// \code
/// VectorValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 2.0},            // composition values (#/unit cell)
///     {"O", 1./6.},
///     {"Va", 5./6.}});
/// \endcode
///
/// Example: Specifying "param_composition"
/// \code
/// VectorValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 1./6.}});         // composition values (comp_x)
/// \endcode
///
monte::VectorValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::VectorValueMap conditions;
  conditions["temperature"] = to_VectorXd(temperature);
  conditions["mol_composition"] =
      make_mol_composition(composition_converter, comp);
  return conditions;
}

/// \brief Helper for making a conditions VectorValueMap for canonical Monte
///     Carlo calculations, interpreted as an increment
///
/// \param temperature The change in temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate `comp` input and convert to  `comp_n` (dmol_composition).
/// \param comp A map of component names (for change in mol per unit cell
///     composition) or axes names (for change in parametric composition) to
///     value.
///
/// \returns VectorValueMap which contains "temperature" and "mol_composition"
///     (increment).
///
/// Example: Specifying "mol_composition" increment
/// \code
/// VectorValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 0.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", -0.01}});
/// \endcode
///
/// Example: Specifying "param_composition" increment
/// \code
/// VectorValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 0.02}});          // composition values (comp_x)
/// \endcode
///
monte::VectorValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::VectorValueMap conditions;
  conditions["temperature"] = to_VectorXd(temperature);
  conditions["mol_composition"] =
      make_mol_composition_increment(composition_converter, comp);
  return conditions;
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
