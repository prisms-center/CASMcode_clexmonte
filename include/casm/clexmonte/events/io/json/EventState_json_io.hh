#ifndef CASM_clexmonte_events_EventState_json_io
#define CASM_clexmonte_events_EventState_json_io

namespace CASM {
class jsonParser;

namespace clexmonte {
struct EventState;
struct EventData;
struct PrimEventData;

jsonParser &to_json(EventData const &event_data, jsonParser &json,
                    PrimEventData const &prim_event_data);

jsonParser &to_json(PrimEventData const &prim_event_data, jsonParser &json);

jsonParser &to_json(EventState const &event_state, jsonParser &json);

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    PrimEventData const &prim_event_data);

jsonParser &to_json(EventState const &event_state, jsonParser &json,
                    EventData const &event_data,
                    PrimEventData const &prim_event_data);

}  // namespace clexmonte
}  // namespace CASM

#endif
