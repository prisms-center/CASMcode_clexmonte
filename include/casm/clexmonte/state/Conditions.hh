#ifndef CASM_clexmonte_state_Conditions
#define CASM_clexmonte_state_Conditions

#include "casm/clexmonte/state/CorrMatchingPotential.hh"
#include "casm/monte/state/ValueMap.hh"

namespace CASM {

namespace composition {
class CompositionConverter;
}

namespace clexmonte {

/// \brief Holds conditions in form preferrable to monte::ValueMap for
/// calculation
///
/// Notes:
/// - Can also be used to specify a conditions increment when specifying a path
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
  /// Eigen::VectorXd x = (comp_x - *m_condition.param_comp_quad_pot_target());
  /// Eigen::VectorXd const &V = *m_condition.param_comp_quad_pot_matrix();
  /// potential_energy += x.dot(V * x);
  /// \endcode
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
};

/// \brief Return initial + n_increment*increment
Conditions make_incremented_conditions(
    Conditions initial, Conditions const &increment, double n_increment,
    xtal::BasicStructure const &prim,
    composition::CompositionConverter const &composition_converter);

/// \brief Make a monte::ValueMap from Conditions
monte::ValueMap make_value_map_from_conditions(Conditions const &conditions);

/// \brief Make a monte::ValueMap from Conditions representing an increment size
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
