#include "casm/clexmonte/canonical/canonical_impl.hh"
#include "casm/clexmonte/state/make_conditions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

CanonicalPotential::CanonicalPotential(std::shared_ptr<system_type> _system)
    : m_system(_system) {
  if (m_system == nullptr) {
    throw std::runtime_error(
        "Error constructing CanonicalPotential: system is empty");
  }
}

/// \brief Reset pointer to state currently being calculated
///
/// Notes:
/// - If state supercell is modified this must be called again
/// - State DoF values can be modified without calling this again
/// - State conditions can be modified without calling this again
void CanonicalPotential::set(state_type const *state,
                             std::shared_ptr<Conditions> conditions) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting CanonicalPotential state: state is empty");
  }
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

  // conditions-specific
  m_conditions = conditions;
}

/// \brief Pointer to current state
state_type const *CanonicalPotential::state() const { return m_state; }

/// \brief Pointer to current conditions, valid for current state
std::shared_ptr<Conditions> const &CanonicalPotential::conditions() const {
  return m_conditions;
}

/// \brief Pointer to formation energy cluster expansion calculator, valid for
///     current state
std::shared_ptr<clexulator::ClusterExpansion> const &
CanonicalPotential::formation_energy() const {
  return m_formation_energy_clex;
}

/// \brief Calculate (per_supercell) canonical potential value
double CanonicalPotential::per_supercell() {
  return m_formation_energy_clex->per_supercell();
}

/// \brief Calculate (per_unitcell) canonical potential value
double CanonicalPotential::per_unitcell() {
  return m_formation_energy_clex->per_unitcell();
}

/// \brief Calculate change in (per_supercell) canonical potential value due
///     to a series of occupation changes
double CanonicalPotential::occ_delta_per_supercell(
    std::vector<Index> const &linear_site_index,
    std::vector<int> const &new_occ) {
  return m_formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
}

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
///
/// \param temperature The temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert between species per unit cell
///     (mol_composition) and parameteric composition (param_composition).
/// \param comp A map of component names (for species per unit cell composition)
/// or
///     axes names (for parametric composition) to value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "mol_composition".
///
/// Example: Specifying "mol_composition"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 2.0},            // composition values (#/unit cell)
///     {"O", 1./6.},
///     {"Va", 5./6.}});
/// \endcode
///
/// Example: Specifying "param_composition"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 1./6.}});         // composition values (param_composition)
/// \endcode
///
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["mol_composition"] =
      make_mol_composition(composition_converter, comp);
  return conditions;
}

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations, interpreted as an increment
///
/// \param temperature The change in temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert between species per unit cell
///     (mol_composition) and parameteric composition (param_composition).
/// \param comp A map of component names (for change in mol per unit cell
///     composition) or axes names (for change in parametric composition) to
///     value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "mol_composition" (increment).
///
/// Example: Specifying "mol_composition" increment
/// \code
/// ValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 0.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", -0.01}});
/// \endcode
///
/// Example: Specifying "param_composition" increment
/// \code
/// ValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 0.02}});          // composition values (param_composition)
/// \endcode
///
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["mol_composition"] =
      make_mol_composition_increment(composition_converter, comp);
  return conditions;
}

template struct Canonical<std::mt19937_64>;

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
