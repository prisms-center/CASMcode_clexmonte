#ifndef CASM_clexmonte_events_CompleteEventCalculator
#define CASM_clexmonte_events_CompleteEventCalculator

#include "casm/clexmonte/events/EventStateCalculator.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief CompleteEventCalculator is an event calculator with the required
/// interface for the
///     classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct CompleteEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Prim event calculators - order must match prim_event_list
  std::vector<EventStateCalculator> const &prim_event_calculators;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Write to warn about non-normal events
  Log &event_log;

  // Note: to keep all event state calculations, comment out this:
  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Count not-normal events
  Index not_normal_count;

  CompleteEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<EventStateCalculator> const &_prim_event_calculators,
      std::map<EventID, EventData> const &_event_list,
      Log &_event_log = CASM::err_log())
      : prim_event_list(_prim_event_list),
        prim_event_calculators(_prim_event_calculators),
        event_list(_event_list),
        event_log(_event_log),
        not_normal_count(0) {}

  /// \brief Get CASM::monte::OccEvent corresponding to given event ID
  double calculate_rate(EventID const &id) {
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
      ++not_normal_count;
    }

    return event_state.rate;
  }
};

}  // namespace clexmonte
}  // namespace CASM

#endif
