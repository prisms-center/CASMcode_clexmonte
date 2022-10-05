#include "casm/clexmonte/semi_grand_canonical_impl.hh"
#include "casm/clexmonte/state/make_conditions.hh"
#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/results/ResultsAnalysisFunction.hh"
#include "casm/monte/state/ValueMap.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

SemiGrandCanonicalPotential::SemiGrandCanonicalPotential(
    std::shared_ptr<system_type> _system)
    : m_system(_system) {
  if (m_system == nullptr) {
    throw std::runtime_error(
        "Error constructing SemiGrandCanonicalPotential: system is "
        "empty");
  }
}

/// \brief Reset pointer to state currently being calculated
///
/// Notes:
/// - If state supercell is modified this must be called again
/// - State DoF values can be modified without calling this again
/// - If state conditions are modified this must be called again
void SemiGrandCanonicalPotential::set(state_type const *state,
                                      std::shared_ptr<Conditions> conditions) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting SemiGrandCanonicalPotential state: state is empty");
  }
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");
  m_convert = &get_index_conversions(*m_system, *m_state);
  m_n_unitcells = get_transformation_matrix_to_super(*m_state).determinant();

  // conditions-specific
  m_conditions = conditions;
  if (!m_conditions->param_chem_pot.has_value()) {
    throw std::runtime_error(
        "Error setting SemiGrandCanonicalPotential state: no param_chem_pot");
  }
  if (!m_conditions->exchange_chem_pot.has_value()) {
    throw std::runtime_error(
        "Error setting SemiGrandCanonicalPotential state: no "
        "exchange_chem_pot");
  }
}

/// \brief Pointer to current state
state_type const *SemiGrandCanonicalPotential::state() const { return m_state; }

/// \brief Pointer to current conditions
std::shared_ptr<Conditions> const &SemiGrandCanonicalPotential::conditions()
    const {
  return m_conditions;
}

/// \brief Calculate (extensive) semi-grand potential value
double SemiGrandCanonicalPotential::extensive_value() {
  Eigen::VectorXi const &occupation = get_occupation(*m_state);
  Eigen::VectorXd mol_composition =
      get_composition_calculator(*m_system).mean_num_each_component(occupation);
  Eigen::VectorXd param_composition =
      get_composition_converter(*m_system).param_composition(mol_composition);
  Eigen::VectorXd const &param_chem_pot = *m_conditions->param_chem_pot;

  double formation_energy = m_formation_energy_clex->extensive_value();

  double potential_energy =
      formation_energy - m_n_unitcells * param_chem_pot.dot(param_composition);

  return potential_energy;
}

/// \brief Calculate change in (extensive) semi-grand potential value due
///     to a series of occupation changes
double SemiGrandCanonicalPotential::occ_delta_extensive_value(
    std::vector<Index> const &linear_site_index,
    std::vector<int> const &new_occ) {
  Eigen::MatrixXd const &exchange_chem_pot = *m_conditions->exchange_chem_pot;
  monte::Conversions const &convert = *m_convert;
  Eigen::VectorXi const &occupation = get_occupation(*m_state);

  double delta_formation_energy =
      m_formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
  double delta_potential_energy = delta_formation_energy;
  for (Index i = 0; i < linear_site_index.size(); ++i) {
    Index l = linear_site_index[i];
    Index asym = convert.l_to_asym(l);
    Index curr_species = convert.species_index(asym, occupation(l));
    Index new_species = convert.species_index(asym, new_occ[i]);
    delta_potential_energy -= exchange_chem_pot(new_species, curr_species);
  }

  return delta_potential_energy;
}

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

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM
