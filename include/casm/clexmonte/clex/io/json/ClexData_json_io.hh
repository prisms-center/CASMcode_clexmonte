#ifndef CASM_clexmonte_clex_ClexData_json_io
#define CASM_clexmonte_clex_ClexData_json_io

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {
struct ClexData;

void parse(InputParser<clexmonte::ClexData> &parser);

}  // namespace clexmonte
}  // namespace CASM

#endif
