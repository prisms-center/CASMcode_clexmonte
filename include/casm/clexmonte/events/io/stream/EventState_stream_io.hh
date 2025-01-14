#ifndef CASM_clexmonte_events_EventState_stream_io
#define CASM_clexmonte_events_EventState_stream_io

#include <iostream>

namespace CASM {

class Log;

namespace clexmonte {

struct EventState;
struct EventData;
struct PrimEventData;

void print(Log &log, EventState const &event_state);

void print(Log &log, EventState const &event_state,
           PrimEventData const &prim_event_data);

void print(Log &log, EventState const &event_state, EventData const &event_data,
           PrimEventData const &prim_event_data);

}  // namespace clexmonte
}  // namespace CASM

#endif
