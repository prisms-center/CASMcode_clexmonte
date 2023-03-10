#ifndef CASM_clexmonte_kinetic_events
#define CASM_clexmonte_kinetic_events

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical_potential.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

/// \brief Data calculated for a single event in a single state
struct EventState {
  bool is_allowed;  ///< Is allowed given current configuration
  double dE_final;  ///< Final state energy, relative to initial state
  double rate;      ///< Occurance rate
};

/// \brief CompleteEventCalculator is an event calculator with the required
///     interface for the classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct CompleteEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Potential
  std::shared_ptr<SemiGrandCanonicalPotential> potential;

  CompleteEventCalculator(
      std::shared_ptr<SemiGrandCanonicalPotential> _potential,
      std::vector<PrimEventData> const &_prim_event_list,
      std::map<EventID, EventData> const &_event_list);

  /// \brief Get CASM::monte::OccEvent corresponding to given event ID
  double calculate_rate(EventID const &id);
};

struct SemiGrandCanonicalEventData {
  SemiGrandCanonicalEventData(
      std::shared_ptr<system_type> system, state_type const &state,
      monte::OccLocation const &occ_location,
      std::vector<monte::OccSwap> const &grand_canonical_swaps,
      std::shared_ptr<SemiGrandCanonicalPotential> potential);

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Calculator for KMC event selection
  std::shared_ptr<CompleteEventCalculator> event_calculator;
};

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
