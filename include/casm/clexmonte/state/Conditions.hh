#ifndef CASM_clexmonte_state_Conditions
#define CASM_clexmonte_state_Conditions

#include "casm/clexmonte/state/CorrMatchingPotential.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/monte/ValueMap.hh"

namespace CASM {

namespace composition {
class CompositionConverter;
}

namespace clexmonte {

struct TemperatureConditionsMixin {
  TemperatureConditionsMixin(double _temperature_is_zero_tol = 1e-10)
      : temperature_is_zero_tol(_temperature_is_zero_tol) {
    set_temperature(0.0);
  }

  /// Tolerance for checking if temperature is zero
  double temperature_is_zero_tol;

  /// \brief Temperature (K)
  double temperature;

  /// \brief 1/(CASM::KB*temperature)
  ///
  /// Note: Use set_temperature to be consistent
  double beta;

  /// \brief Set temperature and beta consistently
  void set_temperature(double _temperature) {
    this->temperature = _temperature;
    if (almost_zero(this->temperature, this->temperature_is_zero_tol)) {
      this->beta = std::numeric_limits<double>::infinity();
    } else {
      this->beta = 1.0 / (CASM::KB * this->temperature);
    }
  }

  void set_temperature(monte::ValueMap const &map) {
    this->set_temperature(map.scalar_values.at("temperature"));
  }

  void put_temperature(monte::ValueMap &map) const {
    map.scalar_values.emplace("temperature", this->temperature);
  }
};

struct CompositionConditionsMixin {
  CompositionConditionsMixin(
      composition::CompositionConverter _composition_converter, double _tol)
      : composition_converter(_composition_converter), tol(_tol) {
    bool is_increment = false;
    this->set_param_composition(
        Eigen::VectorXd::Zero(composition_converter.independent_compositions()),
        is_increment);
  }

  /// Tolerance for comparing compositions
  double tol;

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  /// Parametric composition (depends on composition axes definition)
  Eigen::VectorXd param_composition;

  /// Mol composition, in number per primitive cell
  ///
  /// Note: Use set_param_composition / set_mol_composition to be consistent
  Eigen::VectorXd mol_composition;

  /// \brief Set param_composition and mol_composition consistently
  void set_param_composition(Eigen::VectorXd const &_param_composition,
                             bool is_increment) {
    if (_param_composition.size() !=
        this->composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in CompositionConditionsMixin::set_param_composition: "
          "dimensions mismatch");
    }
    if (!is_increment) {
      this->param_composition = _param_composition;
      this->mol_composition =
          this->composition_converter.mol_composition(this->param_composition);
    } else {
      this->param_composition = _param_composition;
      this->mol_composition =
          this->composition_converter.dmol_composition(this->param_composition);
    }
  }

  /// \brief Set mol_composition and param_composition consistently
  void set_mol_composition(Eigen::VectorXd const &_mol_composition,
                           bool is_increment) {
    if (_mol_composition.size() !=
        this->composition_converter.components().size()) {
      throw std::runtime_error(
          "Error in CompositionConditionsMixin::set_mol_composition: "
          "dimensions mismatch");
    }
    if (!is_increment) {
      this->mol_composition = _mol_composition;
      this->param_composition =
          this->composition_converter.param_composition(this->mol_composition);
    } else {
      this->mol_composition = _mol_composition;
      this->param_composition =
          this->composition_converter.dparam_composition(this->mol_composition);
    }
  }

  /// \brief Set mol_composition and param_composition and check for consistency
  void set_composition(std::optional<Eigen::VectorXd> _param_composition,
                       std::optional<Eigen::VectorXd> _mol_composition,
                       bool is_increment) {
    if (_param_composition.has_value()) {
      if (_param_composition->size() !=
          this->composition_converter.independent_compositions()) {
        throw std::runtime_error(
            "Error in CompositionConditionsMixin::set_composition: "
            "param_composition dimensions mismatch");
      }
    }
    if (_mol_composition.has_value()) {
      if (_mol_composition->size() !=
          this->composition_converter.components().size()) {
        throw std::runtime_error(
            "Error in CompositionConditionsMixin::set_composition: "
            "mol_composition dimensions mismatch");
      }
    }

    if (_mol_composition.has_value()) {
      if (_param_composition.has_value()) {
        // both mol_composition and param_composition
        Eigen::VectorXd equivalent_param_composition;
        if (!is_increment) {
          equivalent_param_composition =
              this->composition_converter.param_composition(*_mol_composition);
        } else {
          equivalent_param_composition =
              this->composition_converter.dparam_composition(*_mol_composition);
        }
        if (!CASM::almost_equal(equivalent_param_composition,
                                *_param_composition, tol)) {
          std::stringstream msg;
          msg << "Error in CompositionConditionsMixin::get_composition: "
              << "'mol_composition' and 'param_composition' are not equivalent."
              << std::endl;
          msg << "mol_composition: " << *_mol_composition << std::endl;
          msg << "equivalent param_composition: "
              << equivalent_param_composition << std::endl;
          msg << "param_composition: " << *_param_composition << std::endl;
          msg << "is_increment: " << is_increment << std::endl;
          throw std::runtime_error(msg.str());
        }
        this->set_param_composition(*_param_composition, is_increment);
        return;
      } else {
        // only mol_composition
        this->set_mol_composition(*_mol_composition, is_increment);
        return;
      }
    } else if (!_param_composition.has_value()) {
      // neither mol_composition nor param_composition
      throw std::runtime_error(
          "Error in CompositionConditionsMixin::get_composition: "
          "one of param_composition or mol_composition must be provided");
    } else {
      // only param_composition
      this->set_param_composition(*_param_composition, is_increment);
      return;
    }
  }

  void set_composition(monte::ValueMap const &map, bool is_increment) {
    std::optional<Eigen::VectorXd> _param_composition;
    std::optional<Eigen::VectorXd> _mol_composition;
    if (map.vector_values.count("param_composition")) {
      _param_composition = map.vector_values.at("param_composition");
    }
    if (map.vector_values.count("mol_composition")) {
      _mol_composition = map.vector_values.at("mol_composition");
    }
    set_composition(_param_composition, _mol_composition, is_increment);
  }

  void put_composition(monte::ValueMap &map) const {
    map.vector_values.emplace("param_composition", this->param_composition);
    map.vector_values.emplace("mol_composition", this->mol_composition);
  }
};

struct OptionalCompositionConditionsMixin {
  OptionalCompositionConditionsMixin(
      composition::CompositionConverter _composition_converter, double _tol)
      : composition_converter(_composition_converter), tol(_tol) {}

  /// Tolerance for comparing compositions
  double tol;

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  /// Parametric composition (depends on composition axes definition)
  std::optional<Eigen::VectorXd> param_composition;

  /// Mol composition, in number per primitive cell
  ///
  /// Note: Use set_param_composition / set_mol_composition to be consistent
  std::optional<Eigen::VectorXd> mol_composition;

  /// \brief Set param_composition and mol_composition consistently
  void set_param_composition(std::optional<Eigen::VectorXd> _param_composition,
                             bool is_increment) {
    if (!_param_composition.has_value()) {
      this->reset_composition();
      return;
    }
    if (_param_composition->size() !=
        this->composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in "
          "OptionalCompositionConditionsMixin::set_param_composition: "
          "dimensions mismatch");
    }
    if (!is_increment) {
      this->param_composition = _param_composition;
      this->mol_composition =
          this->composition_converter.mol_composition(*this->param_composition);
    } else {
      this->param_composition = _param_composition;
      this->mol_composition = this->composition_converter.dmol_composition(
          *this->param_composition);
    }
  }

  /// \brief Set mol_composition and param_composition consistently
  void set_mol_composition(std::optional<Eigen::VectorXd> _mol_composition,
                           bool is_increment) {
    if (!_mol_composition.has_value()) {
      this->reset_composition();
      return;
    }
    if (_mol_composition->size() !=
        this->composition_converter.components().size()) {
      throw std::runtime_error(
          "Error in OptionalCompositionConditionsMixin::set_mol_composition: "
          "dimensions mismatch");
    }
    if (!is_increment) {
      this->mol_composition = _mol_composition;
      this->param_composition =
          this->composition_converter.param_composition(*this->mol_composition);
    } else {
      this->mol_composition = _mol_composition;
      this->param_composition = this->composition_converter.dparam_composition(
          *this->mol_composition);
    }
  }

  /// \brief Set mol_composition and param_composition and check for consistency
  void set_composition(std::optional<Eigen::VectorXd> _param_composition,
                       std::optional<Eigen::VectorXd> _mol_composition,
                       bool is_increment) {
    if (_param_composition.has_value()) {
      if (_param_composition->size() !=
          this->composition_converter.independent_compositions()) {
        throw std::runtime_error(
            "Error in OptionalCompositionConditionsMixin::set_composition: "
            "param_composition dimensions mismatch");
      }
    }
    if (_mol_composition.has_value()) {
      if (_mol_composition->size() !=
          this->composition_converter.components().size()) {
        throw std::runtime_error(
            "Error in OptionalCompositionConditionsMixin::set_composition: "
            "mol_composition dimensions mismatch");
      }
    }

    if (_mol_composition.has_value()) {
      if (_param_composition.has_value()) {
        // both mol_composition and param_composition
        Eigen::VectorXd equivalent_param_composition;
        if (!is_increment) {
          equivalent_param_composition =
              this->composition_converter.param_composition(*_mol_composition);
        } else {
          equivalent_param_composition =
              this->composition_converter.dparam_composition(*_mol_composition);
        }
        if (!CASM::almost_equal(equivalent_param_composition,
                                *_param_composition, tol)) {
          std::stringstream msg;
          msg << "Error in "
                 "OptionalCompositionConditionsMixin::set_composition: "
              << "'mol_composition' and 'param_composition' are not equivalent."
              << std::endl;
          msg << "mol_composition: " << *_mol_composition << std::endl;
          msg << "equivalent param_composition: "
              << equivalent_param_composition << std::endl;
          msg << "param_composition: " << *_param_composition << std::endl;
          msg << "is_increment: " << is_increment << std::endl;
          throw std::runtime_error(msg.str());
        }
        this->set_param_composition(_param_composition, is_increment);
        return;
      } else {
        // only mol_composition
        this->set_mol_composition(_mol_composition, is_increment);
        return;
      }
    } else if (!_param_composition.has_value()) {
      // neither mol_composition nor param_composition
      this->reset_composition();
      return;
    } else {
      // only param_composition
      this->set_param_composition(_param_composition, is_increment);
      return;
    }
  }

  void set_composition(monte::ValueMap const &map, bool is_increment) {
    std::optional<Eigen::VectorXd> _param_composition;
    std::optional<Eigen::VectorXd> _mol_composition;
    if (map.vector_values.count("param_composition")) {
      _param_composition = map.vector_values.at("param_composition");
    }
    if (map.vector_values.count("mol_composition")) {
      _mol_composition = map.vector_values.at("mol_composition");
    }
    set_composition(_param_composition, _mol_composition, is_increment);
  }

  void reset_composition() {
    param_composition.reset();
    mol_composition.reset();
  }

  void put_composition(monte::ValueMap &map) const {
    if (this->param_composition.has_value()) {
      map.vector_values.emplace("param_composition", *this->param_composition);
    }
    if (this->mol_composition.has_value()) {
      map.vector_values.emplace("mol_composition", *this->mol_composition);
    }
  }
};

struct ParamChemPotConditionsMixin {
  ParamChemPotConditionsMixin(
      composition::CompositionConverter const &_composition_converter)
      : composition_converter(_composition_converter) {
    this->set_param_chem_pot(Eigen::VectorXd::Zero(
        composition_converter.independent_compositions()));
  }

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  /// \brief Parameteric chemical potential (conjugate to param_composition)
  ///
  /// potential_energy -= m_condition.param_chem_pot().dot(comp_x)
  Eigen::VectorXd param_chem_pot;

  /// \brief Matrix(new_species, curr_species) of chem_pot(new_species) -
  ///     chem_pot(curr_species)
  ///
  /// \code
  /// delta_potential_energy -=
  ///     conditions.exchange_chem_pot(new_species, curr_species);
  /// \endcode
  ///
  /// Note: Use set_param_chem_pot to be consistent
  Eigen::MatrixXd exchange_chem_pot;

  /// \brief Set param_chem_pot and exchange_chem_pot consistently
  void set_param_chem_pot(Eigen::VectorXd const &_param_chem_pot) {
    if (_param_chem_pot.size() !=
        this->composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in ParamChemPotConditionsMixin::set_param_chem_pot: "
          "dimensions mismatch");
    }
    this->param_chem_pot = _param_chem_pot;
    this->exchange_chem_pot = make_exchange_chemical_potential(
        this->param_chem_pot, this->composition_converter);
  }

  void set_param_chem_pot(monte::ValueMap const &map) {
    this->set_param_chem_pot(map.vector_values.at("param_chem_pot"));
  }

  void put_param_chem_pot(monte::ValueMap &map) const {
    map.vector_values.emplace("param_chem_pot", this->param_chem_pot);
  }
};

struct OptionalParamChemPotConditionsMixin {
  OptionalParamChemPotConditionsMixin(
      composition::CompositionConverter const &_composition_converter)
      : composition_converter(_composition_converter) {}

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  /// \brief Parameteric chemical potential (conjugate to param_composition)
  ///
  /// potential_energy -= m_condition.param_chem_pot().dot(comp_x)
  std::optional<Eigen::VectorXd> param_chem_pot;

  /// \brief Matrix(new_species, curr_species) of chem_pot(new_species) -
  ///     chem_pot(curr_species)
  ///
  /// \code
  /// delta_potential_energy -=
  ///     conditions.exchange_chem_pot(new_species, curr_species);
  /// \endcode
  ///
  /// Note: Use set_param_chem_pot to be consistent
  std::optional<Eigen::MatrixXd> exchange_chem_pot;

  /// \brief Set param_chem_pot and exchange_chem_pot consistently
  void set_param_chem_pot(std::optional<Eigen::VectorXd> _param_chem_pot) {
    this->param_chem_pot = _param_chem_pot;
    if (!this->param_chem_pot.has_value()) {
      this->exchange_chem_pot.reset();
      return;
    }

    if (this->param_chem_pot->size() !=
        this->composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in OptionalParamChemPotConditionsMixin::set_param_chem_pot: "
          "dimensions mismatch");
    }
    this->exchange_chem_pot = make_exchange_chemical_potential(
        *this->param_chem_pot, this->composition_converter);
  }

  void set_param_chem_pot(monte::ValueMap const &map) {
    std::optional<Eigen::VectorXd> _param_chem_pot;
    if (map.vector_values.count("param_chem_pot")) {
      _param_chem_pot = map.vector_values.at("param_chem_pot");
    }
    this->set_param_chem_pot(_param_chem_pot);
  }

  void reset_param_chem_pot() {
    param_chem_pot.reset();
    exchange_chem_pot.reset();
  }

  void put_param_chem_pot(monte::ValueMap &map) const {
    if (this->param_chem_pot.has_value()) {
      map.vector_values.emplace("param_chem_pot", *this->param_chem_pot);
    }
  }
};

struct OptionalParamCompQuadPotConditionsMixin {
  OptionalParamCompQuadPotConditionsMixin(
      composition::CompositionConverter const &_composition_converter)
      : composition_converter(_composition_converter) {}

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  /// \brief Location of quadratic potential min (vector or matrix)
  std::optional<Eigen::VectorXd> param_comp_quad_pot_target;

  /// \brief Quadratic potential coefficients (diagonal terms only)
  ///
  /// \code
  /// Eigen::VectorXd x = (comp_x - *conditions.param_comp_quad_pot_target);
  /// Eigen::VectorXd const &v = *m_condition.param_comp_quad_pot_vector();
  /// potential_energy += v.dot((x.array() * x.array()).matrix());
  /// \endcode
  std::optional<Eigen::VectorXd> param_comp_quad_pot_vector;

  /// \brief Quadratic potential coefficients (full matrix)
  ///
  /// \code
  /// Eigen::VectorXd x = (comp_x -
  /// *m_condition.param_comp_quad_pot_target()); Eigen::VectorXd const &V =
  /// *m_condition.param_comp_quad_pot_matrix(); potential_energy += x.dot(V *
  /// x); \endcode
  std::optional<Eigen::MatrixXd> param_comp_quad_pot_matrix;

  /// \brief Set parameters
  void set_param_comp_quad_pot(
      std::optional<Eigen::VectorXd> _param_comp_quad_pot_target,
      std::optional<Eigen::VectorXd> _param_comp_quad_pot_vector,
      std::optional<Eigen::MatrixXd> _param_comp_quad_pot_matrix) {
    this->param_comp_quad_pot_target = _param_comp_quad_pot_target;
    this->param_comp_quad_pot_vector = _param_comp_quad_pot_vector;
    this->param_comp_quad_pot_matrix = _param_comp_quad_pot_matrix;

    Index n_axes = this->composition_converter.independent_compositions();
    if (this->param_comp_quad_pot_target.has_value()) {
      if (this->param_comp_quad_pot_target->size() != n_axes) {
        throw std::runtime_error(
            "Error in "
            "OptionalParamCompQuadPotConditionsMixin::set_comp_quad_pot: "
            "param_comp_quad_pot_target dimensions mismatch");
      }
    }
    if (this->param_comp_quad_pot_vector.has_value()) {
      if (this->param_comp_quad_pot_vector->size() != n_axes) {
        throw std::runtime_error(
            "Error in "
            "OptionalParamCompQuadPotConditionsMixin::set_comp_quad_pot: "
            "param_comp_quad_pot_vector dimensions mismatch");
      }
    }
    if (this->param_comp_quad_pot_matrix.has_value()) {
      Index n_rows = this->param_comp_quad_pot_matrix->rows();
      Index n_cols = this->param_comp_quad_pot_matrix->cols();
      if (n_rows != n_axes || n_cols != n_axes) {
        throw std::runtime_error(
            "Error in "
            "OptionalParamCompQuadPotConditionsMixin::set_comp_quad_pot: "
            "param_comp_quad_pot_matrix dimensions mismatch");
      }
    }
  }

  void set_param_comp_quad_pot(monte::ValueMap const &map) {
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_target;
    std::optional<Eigen::VectorXd> _param_comp_quad_pot_vector;
    std::optional<Eigen::MatrixXd> _param_comp_quad_pot_matrix;
    if (map.vector_values.count("param_comp_quad_pot_target")) {
      _param_comp_quad_pot_target =
          map.vector_values.at("param_comp_quad_pot_target");
    }
    if (map.vector_values.count("param_comp_quad_pot_vector")) {
      _param_comp_quad_pot_vector =
          map.vector_values.at("param_comp_quad_pot_vector");
    }
    if (map.vector_values.count("param_comp_quad_pot_matrix")) {
      _param_comp_quad_pot_matrix =
          map.vector_values.at("param_comp_quad_pot_matrix");
    }
    this->set_param_comp_quad_pot(_param_comp_quad_pot_target,
                                  _param_comp_quad_pot_vector,
                                  _param_comp_quad_pot_matrix);
  }

  void reset_param_comp_quad_pot() {
    this->param_comp_quad_pot_target.reset();
    this->param_comp_quad_pot_vector.reset();
    this->param_comp_quad_pot_matrix.reset();
  }

  void put_param_comp_quad_pot(monte::ValueMap &map) const {
    if (this->param_comp_quad_pot_target.has_value()) {
      map.vector_values.emplace("param_comp_quad_pot_target",
                                *this->param_comp_quad_pot_target);
    }
    if (this->param_comp_quad_pot_vector.has_value()) {
      map.vector_values.emplace("param_comp_quad_pot_vector",
                                *this->param_comp_quad_pot_vector);
    }
    if (this->param_comp_quad_pot_matrix.has_value()) {
      map.matrix_values.emplace("param_comp_quad_pot_matrix",
                                *this->param_comp_quad_pot_matrix);
    }
  }

  void put_comp_quad_pot(monte::ValueMap &map) const {
    if (this->param_comp_quad_pot_target.has_value()) {
      map.vector_values.emplace("param_comp_quad_pot_target",
                                *this->param_comp_quad_pot_target);
    }
    if (this->param_comp_quad_pot_vector.has_value()) {
      map.vector_values.emplace("param_comp_quad_pot_vector",
                                *this->param_comp_quad_pot_vector);
    }
    if (this->param_comp_quad_pot_matrix.has_value()) {
      map.vector_values.emplace("param_comp_quad_pot_matrix",
                                *this->param_comp_quad_pot_matrix);
    }
  }
};

struct OptionalOrderPotConditionsMixin {
  OptionalOrderPotConditionsMixin(
      composition::CompositionConverter const &_composition_converter)
      : composition_converter(_composition_converter) {}

  /// Defines the parametric composition
  composition::CompositionConverter composition_converter;

  // --- Linear potential in order parameter ---

  /// \brief Linear order parameter potential coefficients
  std::optional<Eigen::VectorXd> order_parameter_pot;

  // --- Quadratic potential in order parameter ---

  std::optional<Eigen::VectorXd> order_parameter_quad_pot_target;
  std::optional<Eigen::VectorXd> order_parameter_quad_pot_vector;
  std::optional<Eigen::MatrixXd> order_parameter_quad_pot_matrix;

  /// \brief Set parameters
  void set_order_parameter_pot(
      std::optional<Eigen::VectorXd> _order_parameter_pot,
      std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target,
      std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector,
      std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix) {
    this->order_parameter_pot = _order_parameter_pot;
    this->order_parameter_quad_pot_target = _order_parameter_quad_pot_target;
    this->order_parameter_quad_pot_vector = _order_parameter_quad_pot_vector;
    this->order_parameter_quad_pot_matrix = _order_parameter_quad_pot_matrix;

    Index n_axes = this->composition_converter.independent_compositions();
    if (this->order_parameter_pot.has_value()) {
      if (this->order_parameter_pot->size() != n_axes) {
        throw std::runtime_error(
            "Error in OptionalOrderPotConditionsMixin::set_comp_quad_pot: "
            "order_parameter_pot dimensions mismatch");
      }
    }
    if (this->order_parameter_quad_pot_target.has_value()) {
      if (this->order_parameter_quad_pot_target->size() != n_axes) {
        throw std::runtime_error(
            "Error in OptionalOrderPotConditionsMixin::set_comp_quad_pot: "
            "order_parameter_quad_pot_target dimensions mismatch");
      }
    }
    if (this->order_parameter_quad_pot_vector.has_value()) {
      if (this->order_parameter_quad_pot_vector->size() != n_axes) {
        throw std::runtime_error(
            "Error in OptionalOrderPotConditionsMixin::set_comp_quad_pot: "
            "order_parameter_quad_pot_vector dimensions mismatch");
      }
    }
    if (this->order_parameter_quad_pot_matrix.has_value()) {
      Index n_rows = this->order_parameter_quad_pot_matrix->rows();
      Index n_cols = this->order_parameter_quad_pot_matrix->cols();
      if (n_rows != n_axes || n_cols != n_axes) {
        throw std::runtime_error(
            "Error in OptionalOrderPotConditionsMixin::set_comp_quad_pot: "
            "order_parameter_quad_pot_matrix dimensions mismatch");
      }
    }
  }

  void set_comp_quad_pot(monte::ValueMap const &map) {
    std::optional<Eigen::VectorXd> _order_parameter_pot;
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_target;
    std::optional<Eigen::VectorXd> _order_parameter_quad_pot_vector;
    std::optional<Eigen::MatrixXd> _order_parameter_quad_pot_matrix;
    if (map.vector_values.count("order_parameter_pot")) {
      _order_parameter_quad_pot_target =
          map.vector_values.at("order_parameter_quad_pot_target");
    }
    if (map.vector_values.count("order_parameter_quad_pot_target")) {
      _order_parameter_quad_pot_target =
          map.vector_values.at("order_parameter_quad_pot_target");
    }
    if (map.vector_values.count("order_parameter_quad_pot_vector")) {
      _order_parameter_quad_pot_vector =
          map.vector_values.at("order_parameter_quad_pot_vector");
    }
    if (map.vector_values.count("order_parameter_quad_pot_matrix")) {
      _order_parameter_quad_pot_matrix =
          map.vector_values.at("order_parameter_quad_pot_matrix");
    }
    this->set_order_parameter_pot(
        _order_parameter_pot, _order_parameter_quad_pot_target,
        _order_parameter_quad_pot_vector, _order_parameter_quad_pot_matrix);
  }

  void reset_order_parameter_pot() {
    this->order_parameter_pot.reset();
    this->order_parameter_quad_pot_target.reset();
    this->order_parameter_quad_pot_vector.reset();
    this->order_parameter_quad_pot_matrix.reset();
  }

  void put_order_parameter_pot(monte::ValueMap &map) const {
    if (this->order_parameter_pot.has_value()) {
      map.vector_values.emplace("order_parameter_pot",
                                *this->order_parameter_pot);
    }
    if (this->order_parameter_quad_pot_target.has_value()) {
      map.vector_values.emplace("order_parameter_quad_pot_target",
                                *this->order_parameter_quad_pot_target);
    }
    if (this->order_parameter_quad_pot_vector.has_value()) {
      map.vector_values.emplace("order_parameter_quad_pot_vector",
                                *this->order_parameter_quad_pot_vector);
    }
    if (this->order_parameter_quad_pot_matrix.has_value()) {
      map.matrix_values.emplace("order_parameter_quad_pot_matrix",
                                *this->order_parameter_quad_pot_matrix);
    }
  }
};

struct IncludeFormationEnergyConditionsMixin {
  IncludeFormationEnergyConditionsMixin() : include_formation_energy(true) {}

  /// Include formation energy?
  bool include_formation_energy;

  /// \brief Set include_formation_energy
  void set_include_formation_energy(bool _include_formation_energy) {
    this->include_formation_energy = _include_formation_energy;
  }

  /// \brief Set include_formation_energy
  void set_include_formation_energy(monte::ValueMap const &map) {
    this->set_include_formation_energy(
        map.boolean_values.at("include_formation_energy"));
  }

  void put_include_formation_energy(monte::ValueMap &map) const {
    map.boolean_values.emplace("include_formation_energy",
                               this->include_formation_energy);
  }
};

struct CorrMatchingPotConditionsMixin {
  CorrMatchingPotConditionsMixin(
      CorrCalculatorFunction _random_alloy_corr_f,
      std::shared_ptr<xtal::BasicStructure const> _prim, double _tol)
      : random_alloy_corr_f(_random_alloy_corr_f), prim(_prim), tol(_tol) {}

  CorrCalculatorFunction random_alloy_corr_f;
  std::shared_ptr<xtal::BasicStructure const> prim;
  double tol;

  // --- Correlations matching potential ---

  std::optional<CorrMatchingParams> corr_matching_pot;

  // --- Random alloy correlations matching potential ---
  std::optional<RandomAlloyCorrMatchingParams> random_alloy_corr_matching_pot;

  /// \brief Set corr_matching_pot
  void set_corr_matching_pot(
      std::optional<CorrMatchingParams> _corr_matching_pot) {
    this->corr_matching_pot = _corr_matching_pot;
  }

  /// \brief Set random_alloy_corr_matching_pot
  void set_random_alloy_corr_matching_pot(
      std::optional<RandomAlloyCorrMatchingParams>
          _random_alloy_corr_matching_pot) {
    this->random_alloy_corr_matching_pot = _random_alloy_corr_matching_pot;
  }

  /// \brief Set corr_matching_pot
  void set_corr_matching_pot(monte::ValueMap const &map) {
    corr_matching_pot =
        ConditionsConstructor<CorrMatchingParams>::from_VectorXd(
            map.vector_values.at("corr_matching_pot"), this->tol);
  }

  /// \brief Set random_alloy_corr_matching_pot
  void set_random_alloy_corr_matching_pot(monte::ValueMap const &map) {
    this->random_alloy_corr_matching_pot =
        ConditionsConstructor<RandomAlloyCorrMatchingParams>::from_VectorXd(
            map.vector_values.at("random_alloy_corr_matching_pot"), *this->prim,
            this->random_alloy_corr_f, this->tol);
  }

  void put_corr_matching_pot(monte::ValueMap &map, bool is_increment) const {
    if (!is_increment) {
      map.vector_values.emplace("corr_matching_pot",
                                to_VectorXd(*this->corr_matching_pot));
    } else {
      map.vector_values.emplace(
          "corr_matching_pot", to_VectorXd_increment(*this->corr_matching_pot));
    }
  }

  void put_random_alloy_corr_matching_pot(monte::ValueMap &map,
                                          bool is_increment) const {
    if (!is_increment) {
      map.vector_values.emplace(
          "random_alloy_corr_matching_pot",
          to_VectorXd(*this->random_alloy_corr_matching_pot));
    } else {
      map.vector_values.emplace(
          "random_alloy_corr_matching_pot",
          to_VectorXd_increment(*this->random_alloy_corr_matching_pot));
    }
  }
};

/// \brief Holds conditions in form preferable to monte::ValueMap for
/// calculation
///
/// Notes:
/// - Can also be used to specify a conditions increment when specifying a
/// path
///   in parameter space
/// - Can be converted to/from a monte::ValueMap which is more convenient for
///   incrementing, etc.
struct Conditions {
  /// Tolerance for comparison operators == and !=
  double tolerance;

  /// \brief Temperature (K)
  double temperature;

  /// \brief 1/(CASM::KB*temperature)
  ///
  /// Note: Use set_temperature to be consistent
  double beta;

  /// Include formation energy?
  bool include_formation_energy;

  // --- Composition ---

  /// Parameteric composition (depends on composition axes definition)
  std::optional<Eigen::VectorXd> param_composition;

  /// Mol composition, in number per primitive cell
  ///
  /// Note: Use set_param_composition / set_mol_composition to be consistent
  std::optional<Eigen::VectorXd> mol_composition;

  // --- Linear potential in param_composition ---

  /// Parameteric chemical potential (conjugate to param_composition)
  ///
  /// potential_energy -= m_condition.param_chem_pot().dot(comp_x)
  std::optional<Eigen::VectorXd> param_chem_pot;

  /// Matrix(new_species, curr_species) of chem_pot(new_species) -
  /// chem_pot(curr_species)
  ///
  /// \code
  /// delta_potential_energy -= conditions.exchange_chem_pot(new_species,
  /// curr_species); \endcode
  std::optional<Eigen::MatrixXd> exchange_chem_pot;

  // --- Quadratic potential in param_composition ---

  /// \brief Location of quadratic potential min (vector or matrix)
  std::optional<Eigen::VectorXd> param_comp_quad_pot_target;

  /// \brief Quadratic potential coefficients (diagonal terms only)
  ///
  /// \code
  /// Eigen::VectorXd x = (comp_x - *conditions.param_comp_quad_pot_target);
  /// Eigen::VectorXd const &v = *m_condition.param_comp_quad_pot_vector();
  /// potential_energy += v.dot((x.array() * x.array()).matrix());
  /// \endcode
  std::optional<Eigen::VectorXd> param_comp_quad_pot_vector;

  /// \brief Quadratic potential coefficients (full matrix)
  ///
  /// \code
  /// Eigen::VectorXd x = (comp_x -
  /// *m_condition.param_comp_quad_pot_target()); Eigen::VectorXd const &V =
  /// *m_condition.param_comp_quad_pot_matrix(); potential_energy += x.dot(V *
  /// x); \endcode
  std::optional<Eigen::MatrixXd> param_comp_quad_pot_matrix;

  // --- Linear potential in order parameter ---

  /// \brief Linear order parameter potential coefficients
  std::optional<Eigen::VectorXd> order_parameter_pot;

  // --- Quadratic potential in order parameter ---

  std::optional<Eigen::VectorXd> order_parameter_quad_pot_target;
  std::optional<Eigen::VectorXd> order_parameter_quad_pot_vector;
  std::optional<Eigen::MatrixXd> order_parameter_quad_pot_matrix;

  // --- Correlations matching potential ---

  std::optional<CorrMatchingParams> corr_matching_pot;

  // --- Random alloy correlations matching potential ---

  std::optional<RandomAlloyCorrMatchingParams> random_alloy_corr_matching_pot;

  /// \brief Construct default Conditions
  Conditions();

  /// \brief Set temperature and beta consistently
  void set_temperature(double _temperature);

  /// \brief Set param_composition and mol_composition consistently, using
  /// param_composition
  void set_param_composition(
      Eigen::VectorXd const &_param_composition,
      composition::CompositionConverter const &_composition_converter,
      bool is_increment);

  /// \brief Set param_composition and mol_composition consistently, using
  /// mol_composition
  void set_mol_composition(
      Eigen::VectorXd const &_mol_composition,
      composition::CompositionConverter const &_composition_converter,
      bool is_increment);

  monte::ValueMap to_value_map(bool is_increment) const;
};

/// \brief Return initial + n_increment*increment
Conditions make_incremented_conditions(
    Conditions initial, Conditions const &increment, double n_increment,
    xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter);

/// \brief Make a monte::ValueMap from Conditions
monte::ValueMap make_value_map_from_conditions(Conditions const &conditions);

/// \brief Make a monte::ValueMap from Conditions representing an increment
/// size
monte::ValueMap make_value_map_from_conditions_increment(
    Conditions const &conditions_increment);

/// \brief Make Conditions from a monte::ValueMap
Conditions make_conditions_from_value_map(
    monte::ValueMap const &map, xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter,
    CorrCalculatorFunction random_alloy_corr_f, double tol);

/// \brief Make Conditions increment from a monte::ValueMap
Conditions make_conditions_increment_from_value_map(
    monte::ValueMap const &map, xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter,
    CorrCalculatorFunction random_alloy_corr_f, double tol);

}  // namespace clexmonte
}  // namespace CASM

#endif
