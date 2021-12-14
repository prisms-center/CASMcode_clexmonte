#ifndef CASM_clexmonte_canonical_InputData_json_io
#define CASM_clexmonte_canonical_InputData_json_io

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {
namespace canonical {

struct InputData;

void parse(InputParser<InputData> &parser);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
