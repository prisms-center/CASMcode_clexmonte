#include "casm/clexmonte/state/Conditions.hh"

#include "casm/composition/CompositionConverter.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/monte/ValueMap.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct default Conditions
///
/// Notes:
/// - Sets tolerance=CASM::TOL
/// - Sets temperature=0.0
/// - Sets include_formation_energy=true
/// - All optional parameters are set to std::nullopt
Conditions::Conditions()
    : tolerance(CASM::TOL), include_formation_energy(true) {
  this->set_temperature(0.0);
}

/// \brief Set temperature and beta consistently
void Conditions::set_temperature(double _temperature) {
  this->temperature = _temperature;
  if (almost_zero(this->temperature, this->tolerance)) {
    this->beta = std::numeric_limits<double>::infinity();
  } else {
    this->beta = 1.0 / (CASM::KB * this->temperature);
  }
}

/// \brief Set param_composition and mol_composition consistently, using
/// param_composition
///
/// \brief _param_composition, Vector of parametric composition, as defined
///     by composition_converter
/// \brief _composition_converter, Converts between number of species per
///     unit cell and parametric composition
/// \param is_increment, If true, treat _param_composition as a change in
/// parametric
///     composition and store the equivalent change in mol_composition
void Conditions::set_param_composition(
    Eigen::VectorXd const &_param_composition,
    composition::CompositionConverter const &_composition_converter,
    bool is_increment) {
  if (!is_increment) {
    this->param_composition = _param_composition;
    this->mol_composition =
        _composition_converter.mol_composition(*this->param_composition);
  } else {
    this->param_composition = _param_composition;
    this->mol_composition =
        _composition_converter.dmol_composition(*this->param_composition);
  }
}

/// \brief Set param_composition and mol_composition consistently, using
/// mol_composition
///
/// \brief _mol_composition, Vector of species per unit cell, as defined
///     by composition_converter
/// \brief _composition_converter, Converts between number of species per
///     unit cell and parametric composition
/// \param is_increment, If true, treat _mol_composition as a change in
/// species per
///     unit cell and store the equivalent change in param_composition
void Conditions::set_mol_composition(
    Eigen::VectorXd const &_mol_composition,
    composition::CompositionConverter const &_composition_converter,
    bool is_increment) {
  if (!is_increment) {
    this->mol_composition = _mol_composition;
    this->param_composition =
        _composition_converter.param_composition(*this->mol_composition);
  } else {
    this->mol_composition = _mol_composition;
    this->param_composition =
        _composition_converter.dparam_composition(*this->mol_composition);
  }
}

monte::ValueMap Conditions::to_value_map(bool is_increment) const {
  if (!is_increment) {
    return make_value_map_from_conditions(*this);
  } else {
    return make_value_map_from_conditions_increment(*this);
  }
}

/// \brief Return initial + n_increment*increment
///
/// Notes:
/// - This make use of the conversions of conditions and conditions
///   increment to/from monte::ValueMap in order to generate incremented
///   conditions
Conditions make_incremented_conditions(
    Conditions initial, Conditions const &increment, double n_increment,
    xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter) {
  monte::ValueMap initial_values = make_value_map_from_conditions(initial);
  monte::ValueMap increment_values =
      make_value_map_from_conditions_increment(increment);
  if (is_mismatched(initial_values, increment_values)) {
    throw std::runtime_error(
        "Error in clexmonte::make_incremented_conditions: Mismatched initial "
        "conditions and conditions increment");
  }
  monte::ValueMap incremented =
      make_incremented_values(initial_values, increment_values, n_increment);
  CorrCalculatorFunction random_alloy_corr_f;
  if (initial.random_alloy_corr_matching_pot.has_value()) {
    random_alloy_corr_f =
        initial.random_alloy_corr_matching_pot->random_alloy_corr_f;
  }
  double tolerance = initial.tolerance;
  return make_conditions_from_value_map(
      incremented, prim, composition_converter, random_alloy_corr_f, tolerance);
}

namespace {

/// \brief Populate a monte::ValueMap from Conditions
///
/// This handles the portion that is the same for conditions and a conditions
/// increment
void _make_value_map_from_conditions(monte::ValueMap &map,
                                     Conditions const &conditions) {
  auto _put_vector = [&](std::optional<Eigen::VectorXd> const &param,
                         std::string name) {
    if (param.has_value()) {
      map.vector_values.emplace(name, *param);
    }
  };
  auto _put_matrix = [&](std::optional<Eigen::MatrixXd> const &param,
                         std::string name) {
    if (param.has_value()) {
      map.matrix_values.emplace(name, *param);
    }
  };

  Conditions const &c = conditions;
  map.scalar_values.emplace("temperature", c.temperature);
  map.boolean_values.emplace("include_formation_energy",
                             c.include_formation_energy);

  _put_vector(c.param_composition, "param_composition");
  _put_vector(c.mol_composition, "mol_composition");
  _put_vector(c.param_chem_pot, "param_chem_pot");

  _put_vector(c.param_comp_quad_pot_target, "param_comp_quad_pot_target");
  _put_vector(c.param_comp_quad_pot_vector, "param_comp_quad_pot_vector");
  _put_matrix(c.param_comp_quad_pot_matrix, "param_comp_quad_pot_matrix");

  _put_vector(c.order_parameter_pot, "order_parameter_pot");
  _put_vector(c.order_parameter_quad_pot_target,
              "order_parameter_quad_pot_target");
  _put_vector(c.order_parameter_quad_pot_vector,
              "order_parameter_quad_pot_vector");
  _put_matrix(c.order_parameter_quad_pot_matrix,
              "order_parameter_quad_pot_matrix");
}

}  // namespace

/// \brief Make a monte::ValueMap from Conditions
monte::ValueMap make_value_map_from_conditions(Conditions const &conditions) {
  monte::ValueMap map;

  // handles the portion that is the same for conditions and
  // conditions_increment
  _make_value_map_from_conditions(map, conditions);

  Conditions const &c = conditions;

  if (c.corr_matching_pot.has_value()) {
    map.vector_values.emplace("corr_matching_pot",
                              to_VectorXd(*c.corr_matching_pot));
  }

  if (c.random_alloy_corr_matching_pot.has_value()) {
    map.vector_values.emplace("random_alloy_corr_matching_pot",
                              to_VectorXd(*c.random_alloy_corr_matching_pot));
  }
  return map;
}

/// \brief Make a monte::ValueMap from Conditions representing an increment size
monte::ValueMap make_value_map_from_conditions_increment(
    Conditions const &conditions_increment) {
  monte::ValueMap map;

  // handles the portion that is the same for conditions and
  // conditions_increment
  _make_value_map_from_conditions(map, conditions_increment);

  Conditions const &c = conditions_increment;

  if (c.corr_matching_pot.has_value()) {
    map.vector_values.emplace("corr_matching_pot",
                              to_VectorXd_increment(*c.corr_matching_pot));
  }

  if (c.random_alloy_corr_matching_pot.has_value()) {
    map.vector_values.emplace(
        "random_alloy_corr_matching_pot",
        to_VectorXd_increment(*c.random_alloy_corr_matching_pot));
  }
  return map;
}

namespace {
/// \brief Make Conditions from a monte::ValueMap
void _make_conditions_from_value_map(
    Conditions &conditions, monte::ValueMap const &map,
    xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter,
    CorrCalculatorFunction random_alloy_corr_f, double tol) {
  auto _get_vector = [&](std::optional<Eigen::VectorXd> &param,
                         std::string name) {
    if (map.vector_values.count(name)) {
      param = map.vector_values.at(name);
    }
  };
  auto _check_components_sized_vector = [&](std::string name) {
    if (map.vector_values.count(name)) {
      auto const &v = map.vector_values.at(name);
      Index n_comp = composition_converter.components().size();
      if (v.size() != n_comp) {
        std::stringstream ss;
        ss << "Error: conditions \"" << name
           << "\" size does not match the number of system components";
        throw std::runtime_error(ss.str());
      }
    }
  };
  auto _check_comp_axes_sized_vector = [&](std::string name) {
    if (map.vector_values.count(name)) {
      auto const &v = map.vector_values.at(name);
      Index n_axes = composition_converter.independent_compositions();
      if (v.size() != n_axes) {
        std::stringstream ss;
        ss << "Error: conditions \"" << name
           << "\" size does not match the number of independent composition "
              "axes";
        throw std::runtime_error(ss.str());
      }
    }
  };
  auto _get_comp_axes_sized_vector = [&](std::optional<Eigen::VectorXd> &param,
                                         std::string name) {
    _check_comp_axes_sized_vector(name);
    _get_vector(param, name);
  };
  auto _get_matrix = [&](std::optional<Eigen::MatrixXd> &param,
                         std::string name) {
    if (map.matrix_values.count(name)) {
      param = map.matrix_values.at(name);
    }
  };
  auto _check_comp_axes_sized_matrix = [&](std::string name) {
    if (map.matrix_values.count(name)) {
      auto const &M = map.matrix_values.at(name);

      Index n_axes = composition_converter.independent_compositions();
      if (M.rows() != n_axes || M.cols() != n_axes) {
        std::stringstream ss;
        ss << "Error: conditions \"" << name
           << "\" size does not match the number of independent composition "
              "axes";
        throw std::runtime_error(ss.str());
      }
    }
  };
  auto _get_comp_axes_sized_matrix = [&](std::optional<Eigen::MatrixXd> &param,
                                         std::string name) {
    _check_comp_axes_sized_matrix(name);
    _get_matrix(param, name);
  };

  Conditions &c = conditions;
  if (map.scalar_values.count("temperature")) {
    c.set_temperature(map.scalar_values.at("temperature"));
  }
  if (map.boolean_values.count("include_formation_energy")) {
    c.include_formation_energy =
        map.boolean_values.at("include_formation_energy");
  }

  _check_comp_axes_sized_vector("param_composition");
  _check_components_sized_vector("mol_composition");

  _get_comp_axes_sized_vector(c.param_chem_pot, "param_chem_pot");
  if (c.param_chem_pot.has_value()) {
    c.exchange_chem_pot = make_exchange_chemical_potential(
        *c.param_chem_pot, composition_converter);
  }

  _get_comp_axes_sized_vector(c.param_comp_quad_pot_target,
                              "param_comp_quad_pot_target");
  _get_comp_axes_sized_vector(c.param_comp_quad_pot_vector,
                              "param_comp_quad_pot_vector");
  _get_comp_axes_sized_matrix(c.param_comp_quad_pot_matrix,
                              "param_comp_quad_pot_matrix");

  _get_vector(c.order_parameter_pot, "order_parameter_pot");
  _get_vector(c.order_parameter_quad_pot_target,
              "order_parameter_quad_pot_target");
  _get_vector(c.order_parameter_quad_pot_vector,
              "order_parameter_quad_pot_vector");
  _get_matrix(c.order_parameter_quad_pot_matrix,
              "order_parameter_quad_pot_matrix");

  if (map.vector_values.count("corr_matching_pot")) {
    c.corr_matching_pot =
        ConditionsConstructor<CorrMatchingParams>::from_VectorXd(
            map.vector_values.at("corr_matching_pot"), tol);
  }

  if (map.vector_values.count("random_alloy_corr_matching_pot")) {
    if (!random_alloy_corr_f) {
      throw std::runtime_error(
          "Error in clexmonte::make_conditions_from_value_map: random alloy "
          "correlations calculator is requested but empty");
    }
    c.random_alloy_corr_matching_pot =
        ConditionsConstructor<RandomAlloyCorrMatchingParams>::from_VectorXd(
            map.vector_values.at("random_alloy_corr_matching_pot"), prim,
            random_alloy_corr_f, tol);
  }
}
}  // namespace

/// \brief Make Conditions from a monte::ValueMap
Conditions make_conditions_from_value_map(
    monte::ValueMap const &map, xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter,
    CorrCalculatorFunction random_alloy_corr_f, double tol) {
  Conditions c;
  bool is_increment = false;

  _make_conditions_from_value_map(c, map, prim, composition_converter,
                                  random_alloy_corr_f, tol);

  if (map.vector_values.count("param_composition")) {
    c.set_param_composition(map.vector_values.at("param_composition"),
                            composition_converter, is_increment);
    // check consistency
    if (map.vector_values.count("mol_composition") &&
        !CASM::almost_equal(*c.mol_composition,
                            map.vector_values.at("mol_composition"))) {
      throw std::runtime_error(
          "Error in make_conditions_from_value_map: Inconsistent "
          "param_composition and mol_composition values");
    }
  } else if (map.vector_values.count("mol_composition")) {
    c.set_mol_composition(map.vector_values.at("mol_composition"),
                          composition_converter, is_increment);
  }

  return c;
}

/// \brief Make Conditions increment from a monte::ValueMap
Conditions make_conditions_increment_from_value_map(
    monte::ValueMap const &map, xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter,
    CorrCalculatorFunction random_alloy_corr_f, double tol) {
  Conditions c;
  bool is_increment = true;

  _make_conditions_from_value_map(c, map, prim, composition_converter,
                                  random_alloy_corr_f, tol);

  if (map.vector_values.count("param_composition")) {
    c.set_param_composition(map.vector_values.at("param_composition"),
                            composition_converter, is_increment);
    // check consistency
    if (map.vector_values.count("mol_composition") &&
        !CASM::almost_equal(*c.mol_composition,
                            map.vector_values.at("mol_composition"))) {
      throw std::runtime_error(
          "Error in make_conditions_increment_from_value_map: Inconsistent "
          "param_composition and mol_composition values");
    }
  } else if (map.vector_values.count("mol_composition")) {
    c.set_mol_composition(map.vector_values.at("mol_composition"),
                          composition_converter, is_increment);
  }

  return c;
}

}  // namespace clexmonte
}  // namespace CASM
