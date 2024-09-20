#include "casm/clexmonte/events/io/json/PrimEventData_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {

jsonParser &to_json(
    clexmonte::PrimEventData const &data, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>>
        event_system,
    occ_events::OccEventOutputOptions const &options) {
  json.put_obj();
  json["event_type_name"] = data.event_type_name;
  json["equivalent_index"] = data.equivalent_index;
  json["is_forward"] = data.is_forward;
  json["prim_event_index"] = data.prim_event_index;
  to_json(data.event, json["event"], event_system, options);
  json["sites"] = data.sites;
  json["occ_init"] = data.occ_init;
  json["occ_final"] = data.occ_final;
  return json;
}

}  // namespace CASM
