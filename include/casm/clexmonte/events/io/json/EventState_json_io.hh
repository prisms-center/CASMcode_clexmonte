#ifndef CASM_clexmonte_events_EventState_json_io
#define CASM_clexmonte_events_EventState_json_io

namespace CASM {
class jsonParser;
template <typename T>
class InputParser;

namespace clexmonte {
struct EventID;
struct EventState;
struct EventData;
struct PrimEventData;

jsonParser &to_json(clexmonte::EventID const &event_id, jsonParser &json);

void parse(InputParser<clexmonte::EventID> &parser);

void from_json(clexmonte::EventID &event_id, jsonParser const &json);

jsonParser &to_json(EventData const &event_data, jsonParser &json);

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
