#ifndef CASM_clexmonte_events_EventStateCalculator
#define CASM_clexmonte_events_EventStateCalculator

#include <memory>

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace clexmonte {

/// \brief Event rate calculation for a particular KMC event
///
/// EventStateCalculator is used to separate the event calculation from the
/// event definition data in PrimEventData. All symmetrically equivalent
/// events can use the same EventStateCalculator, but a simple approach
/// is to create one for each distinct event associated with the primitive
/// cell.
class EventStateCalculator {
 public:
  // /// \brief Constructor
  // EventStateCalculator(
  //     std::shared_ptr<Conditions> _conditions,
  //     std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex,
  //     std::shared_ptr<clexulator::MultiLocalClusterExpansion> _event_clex,
  //     std::string _name, std::map<std::string, Index> _glossary);

  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name);

  /// \brief Reset pointer to state currently being calculated
  void set(monte::State<Configuration> const *state,
           std::shared_ptr<Conditions> conditions);

  /// \brief Pointer to state currently being calculated
  monte::State<Configuration> const *get() const;

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, EventData const &event_data,
                             PrimEventData const &prim_event_data) const;

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// State to use
  monte::State<Configuration> const *m_state;

  /// Conditions
  std::shared_ptr<Conditions> m_conditions;

  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
  Index m_kra_index;
  Index m_freq_index;
};

/// \brief Construct a vector EventStateCalculator, one per event in a
///     vector of PrimEventData
std::vector<EventStateCalculator> make_prim_event_calculators(
    std::shared_ptr<system_type> system,
    monte::State<clexmonte::Configuration> const &state,
    std::vector<PrimEventData> const &prim_event_list,
    std::shared_ptr<Conditions> conditions);

/// \brief Set potential calculator so it evaluates using `state`
void set(EventStateCalculator &potential,
         monte::State<Configuration> const &state,
         std::shared_ptr<Conditions> conditions);

// --- Implementation ---

// /// \brief Constructor
// ///
// /// \param _conditions Conditions under which the event rate is calculated
// /// \param _formation_energy_clex Formation energy cluster expansion
// /// \param _event_clex Used to calculate event kra and attempt frequency
// /// \param _name Event name (primarily for error messages and debugging)
// /// \param _glossary Map of <key>:<index> giving the index in
// ///     `_event_clex` for the coefficients / resulting values for
// ///     kra (key="kra") and attempt frequency (key="freq") calculations.
// inline EventStateCalculator::EventStateCalculator(
//     std::shared_ptr<Conditions> _conditions,
//     std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex,
//     std::shared_ptr<clexulator::MultiLocalClusterExpansion> _event_clex,
//     std::string _name, std::map<std::string, Index> _glossary)
//     : conditions(_conditions),
//       formation_energy_clex(_formation_energy_clex),
//       event_clex(_event_clex),
//       name(_name) {
//   auto _check_coeffs = [&](Index &coeff_index, std::string key) {
//     if (!_glossary.count(key)) {
//       std::stringstream ss;
//       ss << "Error constructing " << name << " EventStateCalculator: No " <<
//       key
//          << " cluster expansion";
//       throw std::runtime_error(ss.str());
//     }
//     coeff_index = _glossary.at(key);
//     if (coeff_index < 0 || coeff_index >= event_clex->coefficients().size())
//     {
//       std::stringstream ss;
//       ss << "Error constructing " << name << " EventStateCalculator: " << key
//          << " index out of range";
//       throw std::runtime_error(ss.str());
//     }
//   };
//   _check_coeffs(kra_index, "kra");
//   _check_coeffs(freq_index, "freq");
// }

/// \brief Constructor
inline EventStateCalculator::EventStateCalculator(
    std::shared_ptr<system_type> _system, std::string _event_type_name)
    : m_system(_system), m_event_type_name(_event_type_name) {}

/// \brief Reset pointer to state currently being calculated
inline void EventStateCalculator::set(monte::State<Configuration> const *state,
                                      std::shared_ptr<Conditions> conditions) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting EventStateCalculator state: state is empty");
  }
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");
  // m_formation_energy_clex->set(&get_dof_values(*m_state));

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

/// \brief Pointer to state currently being calculated
inline monte::State<Configuration> const *EventStateCalculator::get() const {
  return m_state;
}

/// \brief Calculate the state of an event
///
/// \param state Stores whether the event is allowed, is "normal",
///     energy barriers, and event rate
/// \param conditions Holds any conditions (i.e. beta) necessary for
///     calculating event state
/// \param event_data Holds information about the particular translational
///     instance of the event, such as linear sites indices and linear
///     unitcell index, necessary for calculating event state.
/// \param event_data Holds information about the event that does not
///     depend on the particular translational instance, such as the
///     initial and final occupation variables.
inline void EventStateCalculator::calculate_event_state(
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
inline std::vector<EventStateCalculator> make_prim_event_calculators(
    std::shared_ptr<system_type> system,
    monte::State<clexmonte::Configuration> const &state,
    std::vector<PrimEventData> const &prim_event_list,
    std::shared_ptr<Conditions> conditions) {
  std::vector<EventStateCalculator> prim_event_calculators;
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(system,
                                        prim_event_data.event_type_name);
    set(prim_event_calculators.back(), state, conditions);
  }
  return prim_event_calculators;
}

/// \brief Set potential calculator so it evaluates using `state`
inline void set(EventStateCalculator &prim_event_calculator,
                monte::State<Configuration> const &state,
                std::shared_ptr<Conditions> conditions) {
  prim_event_calculator.set(&state, conditions);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
