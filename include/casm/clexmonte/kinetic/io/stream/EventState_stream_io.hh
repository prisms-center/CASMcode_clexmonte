#ifndef CASM_clexmonte_kinetic_EventState_stream_io
#define CASM_clexmonte_kinetic_EventState_stream_io

#include <iostream>

namespace CASM {
namespace clexmonte {
struct EventData;
struct PrimEventData;

namespace kinetic {
struct EventState;

void print(std::ostream &out, kinetic::EventState const &event_state);

void print(std::ostream &out, kinetic::EventState const &event_state,
           PrimEventData const &prim_event_data);

void print(std::ostream &out, kinetic::EventState const &event_state,
           EventData const &event_data, PrimEventData const &prim_event_data);

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
