#include "casm/clexmonte/canonical/conditions.hh"

#include "casm/composition/CompositionConverter.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

namespace make_conditions_impl {

monte::VectorValueMap _make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp, bool is_increment) {
  monte::VectorValueMap conditions;
  std::string err_msg;
  if (is_increment) {
    err_msg =
        "Error in CASM::clexmonte::canonical::make_conditions_increment: ";
  } else {
    err_msg = "Error in CASM::clexmonte::canonical::make_conditions: ";
  }

  // set temperature
  Eigen::VectorXd vector_temperature(1);
  vector_temperature(0) = temperature;
  conditions["temperature"] = vector_temperature;

  // attempt to read composition
  bool is_comp_n = false;
  bool is_comp_x = false;

  // input may be mol_composition ("comp_n") or param_composition ("comp_x")
  std::vector<std::string> components = composition_converter.components();
  Eigen::VectorXd vector_comp_n = Eigen::VectorXd::Zero(components.size());
  std::vector<std::string> axes = composition_converter.axes();
  Eigen::VectorXd vector_comp_x = Eigen::VectorXd::Zero(axes.size());
  for (auto element : comp) {
    // is key a component name?
    auto it = std::find(components.begin(), components.end(), element.first);
    if (it != components.end()) {
      is_comp_n = true;
      vector_comp_n(std::distance(components.begin(), it)) = element.second;
      continue;
    }

    // is key an axis name?
    it = std::find(axes.begin(), axes.end(), element.first);
    if (it != axes.end()) {
      is_comp_x = true;
      vector_comp_x(std::distance(axes.begin(), it)) = element.second;
      continue;
    }

    // if not found in components or axes, then there is an error
    is_comp_n = false;
    is_comp_x = false;
    break;
  }

  if (is_comp_n == is_comp_x) {
    std::stringstream msg;
    msg << err_msg << "Invalid occupant or axes names";
    throw std::runtime_error(msg.str());
  }

  if (is_comp_n) {
    double expected_comp_n_sum =
        is_increment ? 0.0 : composition_converter.origin().sum();
    // check sum is expected value
    // (depending on context, should be number of sublatticese or zero)
    if (!CASM::almost_equal(vector_comp_n.sum(), expected_comp_n_sum, TOL)) {
      std::stringstream msg;
      msg << err_msg << "sum != " << expected_comp_n_sum;
      throw std::runtime_error(msg.str());
    }
  } else {
    if (!is_increment) {
      vector_comp_n = composition_converter.mol_composition(vector_comp_x);
    } else {
      vector_comp_n = composition_converter.dmol_composition(vector_comp_x);
    }
  }

  conditions["comp_n"] = vector_comp_n;

  return conditions;
}
}  // namespace make_conditions_impl

/// \brief Helper for making a conditions VectorValueMap for canonical Monte
///     Carlo calculations
///
/// \param temperature The temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate `comp` input and convert to  `comp_n` (mol_composition).
/// \param comp A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
///
/// \returns VectorValueMap which contains "temperature" and "comp_n".
///
/// Example: Specifying "comp_n"
/// \code
/// VectorValueMap conditions = make_canonical_conditions(
///    300.0,                   // temperature (K)
///    composition_components,  // composition vector order
///    {{"Zr", 2.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", 1.99}});
/// \endcode
///
/// Example: Specifying "comp_x"
/// \code
/// VectorValueMap conditions = make_canonical_conditions(
///    300.0,                   // temperature (K)
///    composition_components,  // composition vector order
///    {{"a", 0.005}});         // composition values (comp_x)
/// \endcode
///
monte::VectorValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  using namespace make_conditions_impl;
  bool is_increment = false;
  return _make_conditions(temperature, composition_converter, comp,
                          is_increment);
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
/// \returns VectorValueMap which contains "temperature" and "comp_n"
///     (increment).
///
/// Example: Specifying "comp_n" increment
/// \code
/// VectorValueMap conditions = make_canonical_conditions(
///    10.0,                    // temperature (K)
///    composition_components,  // composition vector order
///    {{"Zr", 0.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", -0.01}});
/// \endcode
///
/// Example: Specifying "comp_x" increment
/// \code
/// VectorValueMap conditions = make_canonical_conditions(
///    10.0,                    // temperature (K)
///    composition_components,  // composition vector order
///    {{"a", 0.02}});          // composition values (comp_x)
/// \endcode
///
monte::VectorValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  using namespace make_conditions_impl;
  bool is_increment = true;
  return _make_conditions(temperature, composition_converter, comp,
                          is_increment);
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
