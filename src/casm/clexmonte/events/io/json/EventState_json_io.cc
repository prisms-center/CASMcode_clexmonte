#include "casm/clexmonte/events/io/json/EventState_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace clexmonte {

jsonParser &to_json(EventState const &event_state, jsonParser &json) {
  json["is_allowed"] = event_state.is_allowed;
  if (event_state.is_allowed) {
    json["is_normal"] = event_state.is_normal;
    json["dE_final"] = event_state.dE_final;
    json["Ekra"] = event_state.Ekra;
    json["dE_activated"] = event_state.dE_activated;
    json["freq"] = event_state.freq;
    json["rate"] = event_state.rate;
  }
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    PrimEventData const &prim_event_data) {
  json["prim_event_index"] = prim_event_data.prim_event_index;
  json["event_type_name"] = prim_event_data.event_type_name;
  json["equivalent_index"] = prim_event_data.equivalent_index;
  json["is_forward"] = prim_event_data.is_forward;
  json["occ_init"] = prim_event_data.occ_init;
  json["occ_final"] = prim_event_data.occ_final;

  to_json(event_state, json);
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    EventData const &event_data,
                    PrimEventData const &prim_event_data) {
  json["unitcell_index"] = event_data.unitcell_index;
  json["linear_site_index"] = event_data.event.linear_site_index;
  to_json(event_state, json, prim_event_data);
  return json;
}

}  // namespace clexmonte
}  // namespace CASM
