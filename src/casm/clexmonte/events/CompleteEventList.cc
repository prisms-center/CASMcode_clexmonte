#include "casm/clexmonte/events/CompleteEventList.hh"

#include "casm/clexmonte/events/event_methods.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

CompleteEventList make_complete_event_list(
    std::vector<PrimEventData> const &prim_event_list,
    std::vector<EventImpactInfo> const &prim_impact_info_list,
    monte::OccLocation const &occ_location,
    std::vector<EventFilterGroup> const &event_filters) {
  CompleteEventList event_list;

  if (prim_event_list.size() != prim_impact_info_list.size()) {
    throw std::runtime_error(
        "Error in make_complete_event_list: prim_event_list and "
        "prim_impact_info_list size mismatch");
  }

  auto const &unitcell_index_converter =
      occ_location.convert().unitcell_index_converter();
  Index n_unitcells = unitcell_index_converter.total_sites();

  RelativeEventImpactTable relative_impact_table(prim_impact_info_list,
                                                 unitcell_index_converter);

  for (Index unitcell_index = 0; unitcell_index < n_unitcells;
       ++unitcell_index) {
    EventFilterGroup const *filter = nullptr;
    for (auto const &test_filter : event_filters) {
      if (test_filter.unitcell_index.count(unitcell_index)) {
        filter = &test_filter;
        break;
      }
    }

    for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
         ++prim_event_index) {
      if (filter) {
        if (filter->include_by_default == true &&
            filter->prim_event_index.count(prim_event_index)) {
          continue;
        }
        if (filter->include_by_default == false &&
            !filter->prim_event_index.count(prim_event_index)) {
          continue;
        }
      }

      PrimEventData const &prim_event_data = prim_event_list[prim_event_index];

      // set event_id
      EventID event_id;
      event_id.prim_event_index = prim_event_index;
      event_id.unitcell_index = unitcell_index;

      // set event_data
      EventData event_data;
      event_data.unitcell_index = unitcell_index;
      xtal::UnitCell translation = unitcell_index_converter(unitcell_index);
      set_event(event_data.event, prim_event_data, translation, occ_location);

      // add to lists
      event_list.impact_table.emplace(event_id,
                                      relative_impact_table(event_id));
      event_list.events.emplace(event_id, event_data);
    }
  }
  return event_list;
}

/// \brief Construct a vector of all EventID
std::vector<EventID> make_complete_event_id_list(
    Index n_unitcells, std::vector<PrimEventData> const &prim_event_list) {
  std::vector<EventID> event_id_list;
  EventID event_id;
  for (Index unitcell_index = 0; unitcell_index < n_unitcells;
       ++unitcell_index) {
    for (Index prim_event_index = 0; prim_event_index < prim_event_list.size();
         ++prim_event_index) {
      // set event_id
      event_id.prim_event_index = prim_event_index;
      event_id.unitcell_index = unitcell_index;
      event_id_list.push_back(event_id);
    }
  }
  return event_id_list;
}

}  // namespace clexmonte
}  // namespace CASM
