#include "casm/clexmonte/state/make_conditions.hh"

#include "casm/composition/CompositionConverter.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace clexmonte {

namespace make_conditions_impl {

/// \brief Convert dict input to vector values
///
/// Input may be specified by mol ("n") (i.e. component names) or parametric
/// ("x") (i.e. axis names):
///
/// - Example, mol ("n") input: { "Zr": 2.0, "O": 1.0, "Va": 1.0}
///   - Throw if not all components provided
/// - Example, parametric ("x") input: {"a": 0.1666666, "b": 0.5}
///   - Throw if not all axes provided
struct InputToVectors {
  composition::CompositionConverter const &composition_converter;
  bool is_increment;
  bool as_mol;
  bool do_not_convert;
  std::string err_msg;
  std::string name_n;
  std::string name_x;
  bool found_n;
  bool found_x;
  Eigen::VectorXd vector_n;
  Eigen::VectorXd vector_x;

  /// Constructor, parses input map -> vector_n or vector_x
  InputToVectors(
      composition::CompositionConverter const &_composition_converter,
      std::map<std::string, double> input, bool _is_composition,
      bool _is_increment, bool _as_mol, bool _do_not_convert);

  /// \brief Return mol_composition / mol_composition_increment /
  ///     param_composition / param_composition_increment, based on choice of
  ///     `is_increment`, `as_mol`
  Eigen::VectorXd make_comp();

  /// \brief Return param_chem_pot / param_chem_pot_increment,
  ///     based on choice of `is_increment`
  Eigen::VectorXd make_param_chem_pot();
};

InputToVectors::InputToVectors(
    composition::CompositionConverter const &_composition_converter,
    std::map<std::string, double> input, bool _is_composition,
    bool _is_increment, bool _as_mol, bool _do_not_convert)
    : composition_converter(_composition_converter),
      is_increment(_is_increment),
      as_mol(_as_mol),
      do_not_convert(_do_not_convert) {
  if (_is_composition) {
    name_n = "mol_composition";
    name_x = "param_composition";
  } else {
    name_n = "mol_chem_pot";
    name_x = "param_chem_pot";
  }

  if (as_mol) {
    if (is_increment) {
      err_msg = "Error making " + name_n + " increment: ";
    } else {
      err_msg = "Error making " + name_n + ": ";
    }
  } else {
    if (is_increment) {
      err_msg = "Error making " + name_x + " increment: ";
    } else {
      err_msg = "Error making " + name_x + ": ";
    }
  }

  // attempt to read composition
  found_n = false;
  found_x = false;

  // input may be mol ("n") or parametric ("x")
  std::vector<std::string> components = composition_converter.components();
  vector_n = Eigen::VectorXd::Zero(components.size());
  std::vector<std::string> axes = composition_converter.axes();
  vector_x = Eigen::VectorXd::Zero(axes.size());
  for (auto element : input) {
    // is key a component name?
    auto it = std::find(components.begin(), components.end(), element.first);
    if (it != components.end()) {
      found_n = true;
      vector_n(std::distance(components.begin(), it)) = element.second;
      continue;
    }

    // is key an axis name?
    it = std::find(axes.begin(), axes.end(), element.first);
    if (it != axes.end()) {
      found_x = true;
      vector_x(std::distance(axes.begin(), it)) = element.second;
      continue;
    }

    // if not found in components or axes, then there is an error
    found_n = false;
    found_x = false;
    break;
  }

  if (found_n == found_x) {
    std::stringstream msg;
    msg << err_msg << "Invalid occupant or axes names";
    throw std::runtime_error(msg.str());
  }

  if (do_not_convert) {
    if ((as_mol && found_x) || (!as_mol && found_n)) {
      std::stringstream msg;
      msg << err_msg << "Invalid occupant or axes names";
      throw std::runtime_error(msg.str());
    }
  }

  if (found_x) {
    // check for missing axes
    Index n_axes = composition_converter.independent_compositions();
    for (Index i = 0; i < n_axes; ++i) {
      std::string name = composition_converter.comp_var(i);
      if (!input.count(name)) {
        std::stringstream msg;
        msg << err_msg << "Missing axis '" << name << "'";
        throw std::runtime_error(msg.str());
      }
    }
  }

  if (found_n) {
    // check for missing components
    for (Index i = 0; i < components.size(); ++i) {
      std::string name = components[i];
      if (!input.count(name)) {
        std::stringstream msg;
        msg << err_msg << "Missing component '" << name << "'";
        throw std::runtime_error(msg.str());
      }
    }
  }
}

Eigen::VectorXd InputToVectors::make_comp() {
  if (found_n) {
    // check sum is expected value
    // (depending on context, should be number of sublatticese or zero)
    double expected_mol_comp_sum =
        is_increment ? 0.0 : composition_converter.origin().sum();
    if (!CASM::almost_equal(vector_n.sum(), expected_mol_comp_sum, TOL)) {
      std::stringstream msg;
      msg << err_msg << "sum != " << expected_mol_comp_sum;
      throw std::runtime_error(msg.str());
    }
  }

  // If returning composition as mol_composition (number per unit cell)
  if (as_mol) {
    if (found_n) {
      return vector_n;
    } else {
      if (!is_increment) {
        return composition_converter.mol_composition(vector_x);
      } else {
        return composition_converter.dmol_composition(vector_x);
      }
    }
  }
  // If returning composition as param_composition
  else {
    if (found_x) {
      return vector_x;
    } else {
      if (!is_increment) {
        return composition_converter.param_composition(vector_n);
      } else {
        return composition_converter.dparam_composition(vector_n);
      }
    }
  }
}

Eigen::VectorXd InputToVectors::make_param_chem_pot() {
  // If requesting `as_mol==true` (conjugage to number per unit cell)
  if (as_mol || found_n) {
    std::stringstream msg;
    msg << err_msg
        << "chemical potential must be specified in terms of the parametric "
           "composition axes";
    throw std::runtime_error(msg.str());
  }
  if (found_x) {
    return vector_x;
  } else {
    std::stringstream msg;
    msg << err_msg << "parametric chemical potential not found";
    throw std::runtime_error(msg.str());
  }
}

}  // namespace make_conditions_impl

/// \brief Helper for making a conditions VectorValueMap, mol_composition
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
/// \param do_not_convert If true, throw if parameteric composition is given.
///
/// Example: Specifying mol_composition via number per unit cell
/// \code
/// Eigen::VectorXd mol_composition = make_mol_composition(
///    composition_converter,   // composition vector order
///    {{"Zr", 2.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", 1.99}});
/// \endcode
///
/// Example: Specifying mol_composition via parametric composition
/// \code
/// Eigen::VectorXd mol_composition = make_mol_composition(
///    composition_converter,   // composition vector order
///    {{"a", 0.005}});         // composition values (param_composition)
/// \endcode
///
Eigen::VectorXd make_mol_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input, bool do_not_convert) {
  bool is_composition = true;
  bool is_increment = false;
  bool as_mol = true;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_comp();
}

/// \brief Helper for making a conditions VectorValueMap, mol_composition
///     increment
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
/// \param do_not_convert If true, throw if parameteric composition is given.
Eigen::VectorXd make_mol_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input, bool do_not_convert) {
  bool is_composition = true;
  bool is_increment = true;
  bool as_mol = true;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_comp();
}

/// \brief Helper for making a conditions VectorValueMap, param_composition
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
/// \param do_not_convert If true, throw if mol composition is given.
Eigen::VectorXd make_param_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input, bool do_not_convert) {
  bool is_composition = true;
  bool is_increment = false;
  bool as_mol = false;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_comp();
}

/// \brief Helper for making a conditions VectorValueMap, param_comp increment
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of component names (for mol per unit cell composition) or
///     axes names (for parametric composition) to value.
/// \param do_not_convert If true, throw if mol composition is given.
Eigen::VectorXd make_param_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input, bool do_not_convert) {
  bool is_composition = true;
  bool is_increment = true;
  bool as_mol = false;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_comp();
}

// --- Chemical potential ---

/// \brief Helper for making a conditions VectorValueMap, param_chem_pot
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of axes names (for chemical potential
///     conjugate to parametric composition) to value.
///
Eigen::VectorXd make_param_chem_pot(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input) {
  bool is_composition = false;
  bool is_increment = false;
  bool as_mol = false;
  bool do_not_convert = true;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_param_chem_pot();
}

/// \brief Helper for making a conditions VectorValueMap, param_chem_pot
///     increment
///
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert composition.
/// \param input A map of axes names (for chemical potential
///     conjugate to parametric composition) to value.
///
Eigen::VectorXd make_param_chem_pot_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> input) {
  bool is_composition = false;
  bool is_increment = true;
  bool as_mol = false;
  bool do_not_convert = true;
  make_conditions_impl::InputToVectors f(composition_converter, input,
                                         is_composition, is_increment, as_mol,
                                         do_not_convert);
  return f.make_param_chem_pot();
}

}  // namespace clexmonte
}  // namespace CASM
