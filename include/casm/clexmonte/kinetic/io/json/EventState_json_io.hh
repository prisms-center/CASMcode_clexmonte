#ifndef CASM_clexmonte_kinetic_EventState_json_io
#define CASM_clexmonte_kinetic_EventState_json_io

namespace CASM {
class jsonParser;

namespace clexmonte {
struct EventData;
struct PrimEventData;

namespace kinetic {
struct EventState;

jsonParser &to_json(EventState const &event_state, jsonParser &json);

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    PrimEventData const &prim_event_data);

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    EventData const &event_data,
                    PrimEventData const &prim_event_data);

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
