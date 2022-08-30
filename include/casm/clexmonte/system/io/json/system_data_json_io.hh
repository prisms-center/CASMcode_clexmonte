#ifndef CASM_clexmonte_system_data_json_io
#define CASM_clexmonte_system_data_json_io

#include <map>
#include <memory>
#include <string>

namespace CASM {

template <typename T>
class InputParser;

namespace clexulator {
class Clexulator;
}

namespace config {
struct Prim;
}

namespace clexmonte {
struct BasisSetClusterInfo;
struct EquivalentsInfo;

/// \brief Parse BasisSetClusterInfo from a bspecs.json / eci.json file
void parse(
    InputParser<BasisSetClusterInfo> &parser, config::Prim const &prim,
    std::map<std::string, std::shared_ptr<clexulator::Clexulator>> basis_sets);

/// \brief Parse EquivalentsInfo from JSON
void parse(InputParser<EquivalentsInfo> &parser, config::Prim const &prim);

}  // namespace clexmonte
}  // namespace CASM

#endif
