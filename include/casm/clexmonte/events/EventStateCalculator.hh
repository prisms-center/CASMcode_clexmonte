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
struct EventStateCalculator {
  std::shared_ptr<Conditions> conditions;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> event_clex;
  std::string name;
  Index kra_index;
  Index freq_index;

  /// \brief Constructor
  EventStateCalculator(
      std::shared_ptr<Conditions> _conditions,
      std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex,
      std::shared_ptr<clexulator::MultiLocalClusterExpansion> _event_clex,
      std::string _name, std::map<std::string, Index> _glossary);

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, EventData const &event_data,
                             PrimEventData const &prim_event_data) const;
};

// --- Implementation ---

/// \brief Constructor
///
/// \param _conditions Conditions under which the event rate is calculated
/// \param _formation_energy_clex Formation energy cluster expansion
/// \param _event_clex Used to calculate event kra and attempt frequency
/// \param _name Event name (primarily for error messages and debugging)
/// \param _glossary Map of <key>:<index> giving the index in
///     `_event_clex` for the coefficients / resulting values for
///     kra (key="kra") and attempt frequency (key="freq") calculations.
inline EventStateCalculator::EventStateCalculator(
    std::shared_ptr<Conditions> _conditions,
    std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex,
    std::shared_ptr<clexulator::MultiLocalClusterExpansion> _event_clex,
    std::string _name, std::map<std::string, Index> _glossary)
    : conditions(_conditions),
      formation_energy_clex(_formation_energy_clex),
      event_clex(_event_clex),
      name(_name) {
  auto _check_coeffs = [&](Index &coeff_index, std::string key) {
    if (!_glossary.count(key)) {
      std::stringstream ss;
      ss << "Error constructing " << name << " EventStateCalculator: No " << key
         << " cluster expansion";
      throw std::runtime_error(ss.str());
    }
    coeff_index = _glossary.at(key);
    if (coeff_index < 0 || coeff_index >= event_clex->coefficients().size()) {
      std::stringstream ss;
      ss << "Error constructing " << name << " EventStateCalculator: " << key
         << " index out of range";
      throw std::runtime_error(ss.str());
    }
  };
  _check_coeffs(kra_index, "kra");
  _check_coeffs(freq_index, "freq");
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
  clexulator::ConfigDoFValues const *dof_values = formation_energy_clex->get();

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
  state.dE_final = formation_energy_clex->occ_delta_value(
      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate KRA and attempt frequency
  Eigen::VectorXd const &event_values = event_clex->values(
      event_data.unitcell_index, prim_event_data.equivalent_index);
  state.Ekra = event_values[kra_index];
  state.freq = event_values[freq_index];

  // calculate energy in activated state, check if "normal", calculate rate
  state.dE_activated = state.dE_final * 0.5 + state.Ekra;
  state.is_normal =
      (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
  if (state.dE_activated < state.dE_final) state.dE_activated = state.dE_final;
  if (state.dE_activated < 0.0) state.dE_activated = 0.0;
  state.rate = state.freq * exp(-conditions->beta * state.dE_activated);
}

/// \brief Construct a vector EventStateCalculator, one per event in a
///     vector of PrimEventData
template <typename SystemType>
std::vector<EventStateCalculator> make_prim_event_calculators(
    SystemType &system, monte::State<clexmonte::Configuration> const &state,
    std::vector<PrimEventData> const &prim_event_list,
    std::shared_ptr<Conditions> conditions) {
  std::vector<EventStateCalculator> prim_event_calculators;
  for (auto const &prim_event_data : prim_event_list) {
    LocalMultiClexData event_local_multiclex_data =
        get_local_multiclex_data(system, prim_event_data.event_type_name);
    prim_event_calculators.emplace_back(
        conditions, get_clex(system, state, "formation_energy"),
        get_local_multiclex(system, state, prim_event_data.event_type_name),
        prim_event_data.event_type_name,
        event_local_multiclex_data.coefficients_glossary);
  }
  return prim_event_calculators;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
