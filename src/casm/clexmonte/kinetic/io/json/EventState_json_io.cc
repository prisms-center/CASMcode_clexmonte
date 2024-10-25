#include "casm/clexmonte/kinetic/io/json/EventState_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

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
  to_json(prim_event_data, json);
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

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM
