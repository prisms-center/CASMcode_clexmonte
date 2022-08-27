#ifndef CASM_clexmonte_system_System_json_io
#define CASM_clexmonte_system_System_json_io

#include <map>
#include <vector>

namespace CASM {

template <typename T>
class InputParser;

namespace config {
struct Prim;
}

namespace occ_events {
struct OccEventRep;
struct OccSystem;
}  // namespace occ_events

namespace clexmonte {
struct EquivalentsInfo;
struct OccEventTypeData;
struct System;

/// \brief Parse EquivalentsInfo from JSON
void parse(InputParser<EquivalentsInfo> &parser, config::Prim const &prim);

/// \brief Parse System from JSON
void parse(InputParser<System> &parser);

}  // namespace clexmonte
}  // namespace CASM

#endif
