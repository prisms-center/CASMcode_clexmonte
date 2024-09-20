#include "casm/clexmonte/events/io/json/EventState_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace clexmonte {

jsonParser &to_json(clexmonte::EventID const &event_id, jsonParser &json) {
  json["unitcell_index"] = event_id.unitcell_index;
  json["prim_event_index"] = event_id.prim_event_index;
  return json;
}

void parse(InputParser<clexmonte::EventID> &parser) {
  auto ptr = std::make_unique<clexmonte::EventID>();
  clexmonte::EventID &id = *ptr;
  parser.require(id.unitcell_index, "unitcell_index");
  parser.require(id.prim_event_index, "prim_event_index");
  if (parser.valid()) {
    parser.value = std::move(ptr);
  }
}

void from_json(clexmonte::EventID &event_id, jsonParser const &json) {
  InputParser<clexmonte::EventID> parser{json};
  std::stringstream ss;
  ss << "Error: Invalid clexmonte::EventID object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  event_id = std::move(*parser.value);
}

jsonParser &to_json(EventData const &event_data, jsonParser &json) {
  json["unitcell_index"] = event_data.unitcell_index;
  json["linear_site_index"] = event_data.event.linear_site_index;
  return json;
}

jsonParser &to_json(EventData const &event_data, jsonParser &json,
                    PrimEventData const &prim_event_data) {
  to_json(event_data, json);
  to_json(prim_event_data, json);
  return json;
}

jsonParser &to_json(PrimEventData const &prim_event_data, jsonParser &json) {
  json["prim_event_index"] = prim_event_data.prim_event_index;
  json["event_type_name"] = prim_event_data.event_type_name;
  json["equivalent_index"] = prim_event_data.equivalent_index;
  json["is_forward"] = prim_event_data.is_forward;
  json["occ_init"] = prim_event_data.occ_init;
  json["occ_final"] = prim_event_data.occ_final;
  return json;
}

jsonParser &to_json(EventState const &event_state, jsonParser &json) {
  json["is_allowed"] = event_state.is_allowed;
  if (event_state.is_allowed) {
    if (event_state.formation_energy_delta_corr) {
      to_json(*event_state.formation_energy_delta_corr,
              json["formation_energy_delta_corr"], jsonParser::as_array());
    }
    if (event_state.local_corr) {
      to_json(*event_state.local_corr, json["local_corr"],
              jsonParser::as_array());
    }
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

}  // namespace clexmonte
}  // namespace CASM
