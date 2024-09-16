#ifndef CASM_clexmonte_BaseMonteEventData
#define CASM_clexmonte_BaseMonteEventData

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/System.hh"

namespace CASM {
namespace clexmonte {

/// \brief Base class to provide access to event data for a Monte Carlo
/// simulation
class BaseMonteEventData {
 public:
  BaseMonteEventData() = default;
  virtual ~BaseMonteEventData() = default;

  /// The system
  std::shared_ptr<system_type> system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// Get the formation energy coefficients
  virtual clexulator::SparseCoefficients const &formation_energy_coefficients()
      const = 0;

  /// Get the attempt frequency coefficients for a specific event
  virtual clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const = 0;

  /// Get the KRA coefficients for a specific event
  virtual clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const = 0;

  // -- Event list summary info --

  /// The size of the event list
  virtual Index n_events() const = 0;

  /// Return the current total event rate
  virtual double total_rate() const = 0;

  // -- Event list iteration --

  /// Move internal iterator back to beginning of event list
  virtual void rewind() = 0;

  /// Advance internal iterator by one event
  virtual void advance() = 0;

  /// Check if internal iterator is at the end of the event list
  virtual bool is_end() const = 0;

  /// The event ID for the current state of the internal iterator
  virtual EventID const &event_id() const = 0;

  /// The event data for the current state of the internal iterator
  virtual EventData const &event_data() const = 0;

  // -- Event info (accessed by EventID) --

  /// Return the current rate for a specific event
  virtual double event_rate(EventID const &id) const = 0;

  /// Calculate event state data
  virtual EventState const &event_state(EventID const &id) const = 0;

  /// The events that must be updated if the specified event occurs
  virtual std::vector<EventID> const &impact(EventID const &id) const = 0;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
