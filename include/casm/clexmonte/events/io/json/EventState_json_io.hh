#ifndef CASM_clexmonte_events_EventState_json_io
#define CASM_clexmonte_events_EventState_json_io

namespace CASM {
class jsonParer;

namespace clexmonte {
struct EventState;

jsonParser &to_json(EventState const &event_state, jsonParser &json);

}  // namespace clexmonte
}  // namespace CASM

#endif
