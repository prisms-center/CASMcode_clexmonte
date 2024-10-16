#ifndef CASM_clexmonte_events_AllowedEventList
#define CASM_clexmonte_events_AllowedEventList

#include <map>
#include <vector>

#include "casm/clexmonte/events/ImpactTable.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
class PrimNeighborList;
class SuperNeighborList;
}  // namespace clexulator

namespace monte {
class OccLocation;
}  // namespace monte

namespace clexmonte {

struct AllowedEventData {
  /// \brief Whether the event_id is assigned
  bool is_assigned;

  /// \brief The event ID
  ///
  /// This may be in an invalid state if `is_assigned` is false.
  EventID event_id;
};

/// \brief Data structure for KMC storing only the allowed events
///
/// This is designed to work with `RejectionFreeEventSelector` as follows:
///
/// - `AllowedEventList::events` stores the currently allowed events, plus some
///   extra elements to avoid resizing the event rate tree too often.
/// - When `RejectionFreeEventSelector::select_event` is called, an event is
///   selected as an index into `AllowedEventList::events`.
/// - When `RejectionFreeEventSelector::set_impacted_events` is called, the
///   method `AllowedEventList::make_impact_list` is called which:
///
///   - Updates `AllowedEventList::selected_event_id` to the selected event ID.
///   - Uses `AllowedEventList::relative_impact_table` to determine which
///     events are impacted and update `events` to contain the allowed events in
///     the post- selected event application state.
///   - Updates `AllowedEventList::impact_list` to contain the elements of
///     `events` that must be updated before the next event is selected.
///
/// - When `RejectionFreeEventSelector::select_event` is called and updates
///   rates via `rate_calculator_ptr->calculate_rate(event_id)`, any events
///   that are no longer allowed should be freed via `AllowedEventList::free`.
struct AllowedEventList {
  AllowedEventList(
      std::vector<PrimEventData> const &prim_event_list,
      std::vector<EventImpactInfo> const &prim_impact_info_list,
      clexulator::ConfigDoFValues const &dof_values,
      monte::OccLocation const &occ_location,
      std::shared_ptr<clexulator::PrimNeighborList> prim_nlist,
      std::shared_ptr<clexulator::SuperNeighborList> supercell_nlist);

  /// \brief The relative impact table
  RelativeEventImpactTable relative_impact_table;

  /// \brief A mapping of event IDs to linear indices in the `events` vector to
  /// efficiently generate a list of indices into `events` of impacted events.
  /// This only contains currently allowed events.
  std::map<EventID, Index> event_id_to_event_linear_index;

  /// \brief A vector that specifies the currently allowed events
  ///
  /// This also contains a certain number of events that are not allowed,
  /// typically about 10% of the total vector, to avoid having to resize the
  /// event rate tree too often.
  std::vector<AllowedEventData> events;

  /// \brief The available elements of `events` (where `is_allowed` is false)
  std::set<Index> available_events;

  /// \brief The neighbor indices for the sites that change for each event in
  /// the prim_event_list
  ///
  /// `neighbor_index[prim_event_index][site_index]` is the index in the
  /// neighbor list for the `site_index`-th site in the `prim_event_index`-th
  /// event.
  std::vector<std::vector<int>> neighbor_index;

  /// \brief When an event is selected, this stores selected event ID
  EventID selected_event_id;

  /// \brief If true, the event rate tree must be rebuilt before selecting the
  /// next event
  bool rebuild_event_rate_tree;

  /// \brief A list that gets updated based on the selected event to contain the
  /// elements of `events` that are impacted by the selected event
  std::vector<Index> impact_list;

  /// \brief Returns a list {0, 1, 2, ..., events.size()-1}
  std::vector<Index> event_index_list() const;

  /// \brief Returns a list of indices of events that are impacted by the
  /// selected event
  std::vector<Index> const &make_impact_list(Index selected_event_index);

  /// \brief Free an element of `events` for an event which is no longer
  /// allowed
  void free(EventID event_id);
};

}  // namespace clexmonte
}  // namespace CASM

#endif
