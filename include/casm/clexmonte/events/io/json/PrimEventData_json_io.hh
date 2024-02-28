#ifndef CASM_clexmonte_events_PrimEventData_json_io
#define CASM_clexmonte_events_PrimEventData_json_io

#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"

namespace CASM {

class jsonParser;

namespace clexmonte {
struct PrimEventData;
}  // namespace clexmonte

jsonParser &to_json(clexmonte::PrimEventData const &filter, jsonParser &json,
                    occ_events::OccSystem const &event_system,
                    occ_events::OccEventOutputOptions const &options =
                        occ_events::OccEventOutputOptions());

}  // namespace CASM

#endif
