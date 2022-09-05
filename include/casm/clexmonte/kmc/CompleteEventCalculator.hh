#ifndef CASM_clexmonte_kmc_CompleteEventCalculator
#define CASM_clexmonte_kmc_CompleteEventCalculator

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/kmc/PrimEventCalculator.hh"

namespace CASM {
namespace clexmonte {
namespace kmc {

/// \brief Light-weight EventCalculator holds references to external data
/// structures
///
/// Notes:
/// - Expected to be constructed as shared_ptr
struct CompleteEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Prim event calculators - order must match prim_event_list
  std::vector<PrimEventCalculator> const &prim_event_calculators;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Write to warn about non-normal events
  Log &event_log;

  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Count not-normal events
  Index not_normal_count;

  CompleteEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<PrimEventCalculator> const &_prim_event_calculators,
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

}  // namespace kmc
}  // namespace clexmonte
}  // namespace CASM

#endif
