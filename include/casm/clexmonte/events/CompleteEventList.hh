#ifndef CASM_clexmonte_events_CompleteEventList
#define CASM_clexmonte_events_CompleteEventList

#include <map>
#include <vector>

#include "casm/clexmonte/events/ImpactTable.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace monte {
class OccLocation;
}

namespace clexmonte {

struct CompleteEventList {
  std::map<EventID, std::vector<EventID>> impact_table;

  std::map<EventID, EventData> events;
};

CompleteEventList make_complete_event_list(
    std::vector<PrimEventData> const &prim_event_list,
    std::vector<EventImpactInfo> const &prim_impact_info_list,
    monte::OccLocation const &occ_location,
    std::vector<EventFilterGroup> const &event_filters = {});

std::vector<EventID> make_complete_event_id_list(
    Index n_unitcells, std::vector<PrimEventData> const &prim_event_list);

}  // namespace clexmonte
}  // namespace CASM

#endif
