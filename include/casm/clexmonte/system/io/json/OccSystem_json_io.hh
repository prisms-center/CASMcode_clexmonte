#ifndef CASM_clexmonte_system_OccSystem_json_io
#define CASM_clexmonte_system_OccSystem_json_io

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {
struct OccSystem;

/// \brief Parse clexmonte::EquivalentsInfo from JSON
void parse(
      InputParser<EquivalentsInfo> &parser,
      config::Prim const &prim);

/// \brief Parse "events" from JSON
void parse(
      InputParser<OccEventTypeData> &parser,
      config::Prim const &prim);

/// \brief Parse clexmonte::OccSystem from JSON
void parse(InputParser<OccSystem> &parser);

}  // namespace clexmonte
}  // namespace CASM

#endif
