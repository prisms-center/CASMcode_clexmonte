#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
EventStateCalculator::EventStateCalculator(std::shared_ptr<system_type> _system,
                                           std::string _event_type_name)
    : m_system(_system),
      m_event_type_name(_event_type_name),
      m_custom_event_state_calculation(false),
      m_custom_event_state_calculation_f(nullptr) {}

/// \brief Reset pointer to state currently being calculated
void EventStateCalculator::set(state_type const *state) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting EventStateCalculator state: state is empty");
  }
  m_temperature = &m_state->conditions.scalar_values.at("temperature");
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

  // set and validate event clex
  LocalMultiClexData event_local_multiclex_data =
      get_local_multiclex_data(*m_system, m_event_type_name);
  m_event_clex = get_local_multiclex(*m_system, *m_state, m_event_type_name);
  m_event_values.resize(m_event_clex->coefficients().size());
  std::map<std::string, Index> _glossary =
      event_local_multiclex_data.coefficients_glossary;

  auto _check_coeffs = [&](Index &coeff_index, std::string key) {
    if (!_glossary.count(key)) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: No " << key << " cluster expansion";
      throw std::runtime_error(ss.str());
    }
    coeff_index = _glossary.at(key);
    if (coeff_index < 0 || coeff_index >= m_event_clex->coefficients().size()) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: " << key << " index out of range";
      throw std::runtime_error(ss.str());
    }
  };
  _check_coeffs(m_kra_index, "kra");
  _check_coeffs(m_freq_index, "freq");
}

/// \brief Set custom event state calculation function
void EventStateCalculator::set_custom_event_state_calculation(
    CustomEventStateCalculationFunction f) {
  m_custom_event_state_calculation = true;
  m_custom_event_state_calculation_f = f;
}

/// \brief Clear custom event state calculation function
void EventStateCalculator::clear_custom_event_state_calculation() {
  m_custom_event_state_calculation = false;
  m_custom_event_state_calculation_f = nullptr;
}

/// \brief Calculate the state of an event
///
/// If a custom event state calculation function is set, it is called; otherwise
/// the default event state calculation is used.
///
/// \param state
/// \param unitcell_index
/// \param linear_site_index
/// \param prim_event_data
///
void EventStateCalculator::calculate_event_state(
    EventState &state, Index unitcell_index,
    std::vector<Index> const &linear_site_index,
    PrimEventData const &prim_event_data) const {
  // Initialize event state
  state.formation_energy_delta_corr = nullptr;
  state.local_corr = nullptr;

  // Check if event is allowed based on current occupation
  clexulator::ConfigDoFValues const *dof_values =
      m_formation_energy_clex->get();
  state.is_allowed =
      event_is_allowed(linear_site_index, *dof_values, prim_event_data);
  if (!state.is_allowed) {
    state.rate = 0.0;
    return;
  }

  // Calculate event state
  if (this->m_custom_event_state_calculation) {
    // Set current event details for access by custom event state calculation:
    m_unitcell_index = unitcell_index;
    m_linear_site_index = &linear_site_index;
    m_prim_event_data = &prim_event_data;

    // Call custom event state calculation function
    this->m_custom_event_state_calculation_f(std::ref(state), *this);
    return;
  } else {
    this->_default_event_state_calculation(state, unitcell_index,
                                           linear_site_index, prim_event_data);
  }
}

/// \brief Calculate the state of an event
void EventStateCalculator::_default_event_state_calculation(
    EventState &state, Index unitcell_index,
    std::vector<Index> const &linear_site_index,
    PrimEventData const &prim_event_data) const {
  // calculate change in energy to final state
  //  state.dE_final = m_formation_energy_clex->occ_delta_value(
  //      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate change in energy to final state
  // - and save pointer to delta correlations
  state.formation_energy_delta_corr =
      &m_formation_energy_clex->correlations().occ_delta(
          linear_site_index, prim_event_data.occ_final);
  state.dE_final = m_formation_energy_clex->coefficients() *
                   (*state.formation_energy_delta_corr);

  // calculate KRA and attempt frequency
  // - add save pointer to local correlations
  state.local_corr = &m_event_clex->correlations().local(
      unitcell_index, prim_event_data.equivalent_index);
  for (int i = 0; i < m_event_clex->coefficients().size(); ++i) {
    m_event_values(i) = m_event_clex->coefficients()[i] * (*state.local_corr);
  }
  state.Ekra = m_event_values[m_kra_index];
  state.freq = m_event_values[m_freq_index];

  // calculate energy in activated state, check if "normal", calculate rate
  state.dE_activated = state.dE_final * 0.5 + state.Ekra;
  state.is_normal =
      (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
  if (state.dE_activated < state.dE_final) state.dE_activated = state.dE_final;
  if (state.dE_activated < 0.0) state.dE_activated = 0.0;

  state.rate = state.freq * exp(-this->beta() * state.dE_activated);
}

}  // namespace clexmonte
}  // namespace CASM
