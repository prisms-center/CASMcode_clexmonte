#include "casm/clexmonte/events/AllowedEventList.hh"

#include <numeric>

#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

AllowedEventList::AllowedEventList(
    std::vector<PrimEventData> const &prim_event_list,
    std::vector<EventImpactInfo> const &prim_impact_info_list,
    clexulator::ConfigDoFValues const &dof_values,
    monte::OccLocation const &occ_location,
    std::shared_ptr<clexulator::PrimNeighborList> prim_nlist,
    std::shared_ptr<clexulator::SuperNeighborList> supercell_nlist)
    : relative_impact_table(prim_impact_info_list,
                            occ_location.convert().unitcell_index_converter()) {
  if (prim_event_list.size() != prim_impact_info_list.size()) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: prim_event_list and "
        "prim_impact_info_list size mismatch");
  }
  if (prim_nlist == nullptr) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: prim_nlist is nullptr");
  }
  if (supercell_nlist == nullptr) {
    throw std::runtime_error(
        "Error in AllowedEventList constructor: supercell_nlist is nullptr");
  }

  auto const &unitcell_index_converter =
      occ_location.convert().unitcell_index_converter();
  Index n_unitcells = unitcell_index_converter.total_sites();

  // Construct `neighbor_index`
  this->neighbor_index.reserve(prim_event_list.size());
  for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
       ++prim_event_index) {
    PrimEventData const &prim_event_data = prim_event_list[prim_event_index];

    std::vector<int> _neighbor_index;
    for (auto const &site : prim_event_data.sites) {
      _neighbor_index.push_back(prim_nlist->neighbor_index(site));
    }
    this->neighbor_index.push_back(_neighbor_index);
  }

  // Construct `events` list and consistent `event_id_to_event_linear_index`
  std::vector<Index> linear_site_index;
  for (Index unitcell_index = 0; unitcell_index < n_unitcells;
       ++unitcell_index) {
    for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
         ++prim_event_index) {
      PrimEventData const &prim_event_data = prim_event_list[prim_event_index];

      AllowedEventData allowed_event_data;

      // set event_id
      allowed_event_data.event_id.prim_event_index = prim_event_index;
      allowed_event_data.event_id.unitcell_index = unitcell_index;

      // set linear_site_index
      set_event_linear_site_index(linear_site_index, unitcell_index,
                                  this->neighbor_index[prim_event_index],
                                  *supercell_nlist);

      // set `is_allowed`
      bool is_allowed = true;
      int i = 0;
      for (Index l : linear_site_index) {
        if (dof_values.occupation(l) != prim_event_data.occ_init[i]) {
          is_allowed = false;
          break;
        }
        ++i;
      }

      if (is_allowed) {
        allowed_event_data.is_assigned = true;
        this->event_id_to_event_linear_index.emplace(
            allowed_event_data.event_id, this->events.size());

        this->events.push_back(allowed_event_data);
      }
    }
  }

  // Add elements to `events` to avoid resizing the event rate tree too often
  Index n_extra_events =
      static_cast<Index>(std::ceil(this->events.size() / 9.0));
  for (Index i = 0; i < n_extra_events; ++i) {
    AllowedEventData allowed_event_data;
    allowed_event_data.is_assigned = false;
    this->available_events.insert(events.size());
    this->events.push_back(allowed_event_data);
  }
}

/// \brief Returns a list {0, 1, 2, ..., events.size()-1}
std::vector<Index> AllowedEventList::event_index_list() const {
  std::vector<Index> _list(this->events.size());
  std::iota(_list.begin(), _list.end(), 0);
  return _list;
}

/// \brief Returns a list of indices of events that are impacted by the
/// selected event
///
/// This is designed to work with `RejectionFreeEventSelector`.
/// \param selected_event_index
/// \return
std::vector<Index> const &AllowedEventList::make_impact_list(
    Index selected_event_index) {
  // set `selected_event_id`
  this->selected_event_id = this->events[selected_event_index].event_id;

  // set `impact_list` and update `events`
  impact_list.clear();
  std::vector<EventID> impacted_event_ids =
      this->relative_impact_table(this->selected_event_id);
  for (auto const &event_id : impacted_event_ids) {
    auto it = this->event_id_to_event_linear_index.find(event_id);
    if (it != this->event_id_to_event_linear_index.end()) {
      // impacted event is already assigned
      this->impact_list.push_back(it->second);
    } else {
      // impacted event is not assigned
      if (this->available_events.empty()) {
        // no available events, set flag that the event rate tree must be
        // rebuilt
        this->rebuild_event_rate_tree = true;
        this->impact_list.clear();
        return this->impact_list;
      } else {
        // assign an available event
        auto avail_it = this->available_events.begin();

        AllowedEventData &allowed_event_data = this->events[*avail_it];
        allowed_event_data.is_assigned = true;
        allowed_event_data.event_id = event_id;

        this->event_id_to_event_linear_index.emplace(event_id, *avail_it);
        this->impact_list.push_back(*avail_it);
        this->available_events.erase(avail_it);
      }
    }
  }
  return this->impact_list;
}

/// \brief Free an element of `events` for an event which is no longer
/// allowed
void AllowedEventList::free(EventID event_id) {
  auto it = this->event_id_to_event_linear_index.find(event_id);
  if (it == this->event_id_to_event_linear_index.end()) {
    throw std::runtime_error(
        "Error in AllowedEventList::free: event_id not found");
  }
  AllowedEventData &allowed_event_data = this->events[it->second];
  allowed_event_data.is_assigned = false;
  this->available_events.insert(it->second);
  this->event_id_to_event_linear_index.erase(it);
}

}  // namespace clexmonte
}  // namespace CASM
