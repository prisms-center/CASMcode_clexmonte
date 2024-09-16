#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"

#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/events/event_data.hh"

namespace CASM {
namespace clexmonte {

void print(std::ostream &out, EventState const &event_state) {
  out << "is_allowed: " << std::boolalpha << event_state.is_allowed
      << std::endl;
  if (event_state.is_allowed) {
    out << "dE_activated: " << event_state.dE_activated << std::endl;
    out << "dE_final: " << event_state.dE_final << std::endl;
    out << "is_normal: " << std::boolalpha << event_state.is_normal
        << std::endl;
    out << "Ekra: " << event_state.Ekra << std::endl;
    out << "freq: " << event_state.freq << std::endl;
    out << "rate: " << event_state.rate << std::endl;
  }
}

void print(std::ostream &out, EventState const &event_state,
           PrimEventData const &prim_event_data) {
  out << "prim_event_index: " << prim_event_data.prim_event_index << std::endl;
  out << "event_type_name: " << prim_event_data.event_type_name << std::endl;
  out << "equivalent_index: " << prim_event_data.equivalent_index << std::endl;
  out << "is_forward: " << std::boolalpha << prim_event_data.is_forward
      << std::endl;
  out << "occ_init: " << prim_event_data.occ_init << std::endl;
  out << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(out, event_state);
}

void print(std::ostream &out, EventState const &event_state,
           EventData const &event_data, PrimEventData const &prim_event_data) {
  out << "prim_event_index: " << prim_event_data.prim_event_index << std::endl;
  out << "unitcell_index: " << event_data.unitcell_index << std::endl;
  out << "event_type_name: " << prim_event_data.event_type_name << std::endl;
  out << "equivalent_index: " << prim_event_data.equivalent_index << std::endl;
  out << "is_forward: " << std::boolalpha << prim_event_data.is_forward
      << std::endl;
  out << "linear_site_index: " << event_data.event.linear_site_index
      << std::endl;
  out << "occ_init: " << prim_event_data.occ_init << std::endl;
  out << "occ_final: " << prim_event_data.occ_final << std::endl;
  print(out, event_state);
}

}  // namespace clexmonte
}  // namespace CASM
