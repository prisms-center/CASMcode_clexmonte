#ifndef CASM_clexmonte_MonteEventData
#define CASM_clexmonte_MonteEventData

#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

namespace CASM {
namespace clexmonte {

class MonteEventData {
 public:
  MonteEventData(std::shared_ptr<BaseMonteEventData> _pot,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_data(_pot), m_lib(_lib) {}

  ~MonteEventData() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  /// The system
  std::shared_ptr<system_type> system() const { return m_data->system; }

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> const &prim_event_list() const {
    return m_data->prim_event_list;
  }

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> const &prim_impact_info_list() const {
    return m_data->prim_impact_info_list;
  }

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients() const {
    return m_data->formation_energy_coefficients();
  }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const {
    return m_data->freq_coefficients(prim_event_index);
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const {
    return m_data->kra_coefficients(prim_event_index);
  }

  // -- Event list summary info --

  /// The size of the event list
  Index n_events() const { return m_data->n_events(); }

  /// Return the current total event rate
  double total_rate() const { return m_data->total_rate(); }

  // -- Event list iteration --

  /// Move internal iterator back to beginning of event list
  void rewind() const { return m_data->rewind(); }

  /// Advance internal iterator by one event
  void advance() const { return m_data->advance(); }

  /// Check if internal iterator is at the end of the event list
  bool is_end() const { return m_data->is_end(); }

  /// The event ID for the current state of the internal iterator
  EventID const &event_id() const { return m_data->event_id(); }

  /// The event data for the current state of the internal iterator
  EventData const &event_data() const { return m_data->event_data(); }

  // -- Event info (accessed by EventID) --

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const { return m_data->event_rate(id); }

  /// Calculate event state data
  EventState const &event_state(EventID const &id) const {
    return m_data->event_state(id);
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &impact(EventID const &id) const {
    return m_data->impact(id);
  }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
