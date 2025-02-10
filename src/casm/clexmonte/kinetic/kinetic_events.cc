#include "casm/clexmonte/kinetic/kinetic_events.hh"

#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/kinetic/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/system/System.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Constructor
EventStateCalculator::EventStateCalculator(std::shared_ptr<system_type> _system,
                                           std::string _event_type_name)
    : m_system(_system), m_event_type_name(_event_type_name) {}

/// \brief Reset pointer to state currently being calculated
void EventStateCalculator::set(state_type const *state,
                               std::shared_ptr<Conditions> conditions) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting EventStateCalculator state: state is empty");
  }
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

  // set and validate event clex
  LocalMultiClexData event_local_multiclex_data =
      get_local_multiclex_data(*m_system, m_event_type_name);
  m_event_clex = get_local_multiclex(*m_system, *m_state, m_event_type_name);
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

  // conditions-specific
  m_conditions = conditions;
}

/// \brief Pointer to current state
state_type const *EventStateCalculator::state() const { return m_state; }

/// \brief Pointer to current conditions
std::shared_ptr<Conditions> const &EventStateCalculator::conditions() const {
  return m_conditions;
}

/// \brief Calculate the state of an event
///
/// \param state Stores whether the event is allowed, is "normal",
///     energy barriers, and event rate
/// \param event_data Holds information about the particular translational
///     instance of the event, such as linear sites indices and linear
///     unitcell index, necessary for calculating event state.
/// \param prim_event_data Holds information about the event that does not
///     depend on the particular translational instance, such as the
///     initial and final occupation variables.
void EventStateCalculator::calculate_event_state(
    EventState &state, EventData const &event_data,
    PrimEventData const &prim_event_data) const {
  clexulator::ConfigDoFValues const *dof_values =
      m_formation_energy_clex->get();

  int i = 0;
  for (Index l : event_data.event.linear_site_index) {
    if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
      state.is_allowed = false;
      state.rate = 0.0;
      return;
    }
    ++i;
  }
  state.is_allowed = true;

  // calculate change in energy to final state
  state.dE_final = m_formation_energy_clex->occ_delta_value(
      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate KRA and attempt frequency
  Eigen::VectorXd const &event_values = m_event_clex->values(
      event_data.unitcell_index, prim_event_data.equivalent_index);
  state.Ekra = event_values[m_kra_index];
  state.freq = event_values[m_freq_index];

  // calculate energy in activated state, check if "normal", calculate rate
  state.dE_activated = state.dE_final * 0.5 + state.Ekra;
  state.is_normal =
      (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
  if (state.dE_activated < state.dE_final) state.dE_activated = state.dE_final;
  if (state.dE_activated < 0.0) state.dE_activated = 0.0;
  state.rate = state.freq * exp(-m_conditions->beta * state.dE_activated);
}

/// \brief Construct a vector EventStateCalculator, one per event in a
///     vector of PrimEventData
std::vector<EventStateCalculator> make_prim_event_calculators(
    std::shared_ptr<system_type> system, state_type const &state,
    std::vector<PrimEventData> const &prim_event_list,
    std::shared_ptr<Conditions> conditions) {
  std::vector<EventStateCalculator> prim_event_calculators;
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(system,
                                        prim_event_data.event_type_name);
    prim_event_calculators.back().set(&state, conditions);
  }
  return prim_event_calculators;
}

// CompleteEventCalculator

CompleteEventCalculator::CompleteEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    std::map<EventID, EventData> const &_event_list, Log &_event_log)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_log(_event_log),
      abnormal_count(0) {}

/// \brief Get CASM::monte::OccEvent corresponding to given event ID
double CompleteEventCalculator::calculate_rate(EventID const &id) {
  EventData const &event_data = event_list.at(id);
  PrimEventData const &prim_event_data =
      prim_event_list.at(id.prim_event_index);
  // Note: to keep all event state calculations, uncomment this:
  // EventState &event_state = event_data.event_state;
  prim_event_calculators.at(id.prim_event_index)
      .calculate_event_state(event_state, event_data, prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (event_state.is_allowed && !event_state.is_normal) {
    event_log << "---" << std::endl;
    print(event_log.ostream(), event_state, event_data, prim_event_data);
    event_log << std::endl;
    ++abnormal_count;
  }

  return event_state.rate;
}

KineticEventData::KineticEventData(std::shared_ptr<system_type> _system)
    : system(_system) {
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing KineticEventData: no 'formation_energy' clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});
}

/// \brief Update for given state, conditions, and occupants
///
/// This must be called before first use and when changing supercells,
/// conditions, etc. to build the event list and event calculator
void KineticEventData::update(
    state_type const &state, std::shared_ptr<Conditions> conditions,
    monte::OccLocation const &occ_location,
    std::vector<EventFilterGroup> const &event_filters) {
  // These are constructed/re-constructed so cluster expansions point
  // at the current state
  prim_event_calculators = clexmonte::kinetic::make_prim_event_calculators(
      system, state, prim_event_list, conditions);

  // TODO: rejection-clexmonte option does not require impact table
  event_list = clexmonte::make_complete_event_list(
      prim_event_list, prim_impact_info_list, occ_location, event_filters);

  // Construct CompleteEventCalculator
  event_calculator =
      std::make_shared<clexmonte::kinetic::CompleteEventCalculator>(
          prim_event_list, prim_event_calculators, event_list.events);
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM
