#ifndef CASM_clexmonte_system_System_json_io
#define CASM_clexmonte_system_System_json_io

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "casm/global/filesystem.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexulator {
class Clexulator;
}

namespace config {
struct Prim;
}

namespace occ_events {
struct OccEventRep;
struct OccSystem;
}  // namespace occ_events

namespace clexmonte {
struct BasisSetClusterInfo;
struct EquivalentsInfo;
struct OccEventTypeData;
struct System;

/// \brief Parse System from JSON
void parse(InputParser<System> &parser, std::vector<fs::path> search_path,
           bool verbose = false);

}  // namespace clexmonte
}  // namespace CASM

#endif
