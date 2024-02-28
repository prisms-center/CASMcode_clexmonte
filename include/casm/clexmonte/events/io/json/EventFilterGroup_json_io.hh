#ifndef CASM_clexmonte_events_EventFilterGroup_json_io
#define CASM_clexmonte_events_EventFilterGroup_json_io

namespace CASM {

class jsonParser;
template <typename T>
class InputParser;

namespace clexmonte {

struct EventFilterGroup;

}  // namespace clexmonte

jsonParser &to_json(clexmonte::EventFilterGroup const &filter,
                    jsonParser &json);

void parse(InputParser<clexmonte::EventFilterGroup> &parser);

void from_json(clexmonte::EventFilterGroup &filter, jsonParser const &json);

}  // namespace CASM

#endif
