#ifndef CASM_clexmonte_system_data_json_io
#define CASM_clexmonte_system_data_json_io

#include <map>
#include <memory>
#include <string>

namespace CASM {

template <typename T>
class InputParser;
class jsonParser;

namespace xtal {
class BasicStructure;
}

namespace clexulator {
class Clexulator;
}

namespace config {
struct Prim;
}

namespace clexmonte {
struct BasisSetClusterInfo;
struct LocalBasisSetClusterInfo;
struct EquivalentsInfo;
struct LocalOrbitCompositionCalculatorData;

/// \brief Parse BasisSetClusterInfo from a bspecs.json / eci.json file
void parse(InputParser<BasisSetClusterInfo> &parser, config::Prim const &prim);

/// \brief Output minimal "equivalents info" to JSON
jsonParser &to_json(EquivalentsInfo const &equivalents_info, jsonParser &json,
                    xtal::BasicStructure const &prim);

/// \brief Parse EquivalentsInfo from JSON
void parse(InputParser<EquivalentsInfo> &parser, config::Prim const &prim);

/// \brief Parse LocalBasisSetClusterInfo from a bspecs.json / eci.json file
void parse(InputParser<LocalBasisSetClusterInfo> &parser,
           config::Prim const &prim, EquivalentsInfo const &info);

/// \brief Output LocalOrbitCompositionCalculatorData as JSON
jsonParser &to_json(LocalOrbitCompositionCalculatorData const &data,
                    jsonParser &json);

/// \brief Parse LocalOrbitCompositionCalculatorData from JSON
void parse(InputParser<LocalOrbitCompositionCalculatorData> &parser,
           std::string local_basis_set_name);

}  // namespace clexmonte
}  // namespace CASM

#endif
