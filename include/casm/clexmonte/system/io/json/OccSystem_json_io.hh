#ifndef CASM_clexmonte_system_OccSystem_json_io
#define CASM_clexmonte_system_OccSystem_json_io

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {
struct OccSystem;

/// \brief Parse clexmonte::OccSystem from JSON
void parse(InputParser<OccSystem> &parser);

}  // namespace clexmonte
}  // namespace CASM

#endif
