#include "casm/clexmonte/events/io/json/EventFilterGroup_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"

namespace CASM {

jsonParser &to_json(clexmonte::EventFilterGroup const &filter,
                    jsonParser &json) {
  json.put_obj();
  json["unitcell_index"] = filter.unitcell_index;
  json["include_by_default"] = filter.include_by_default;
  json["prim_event_index"] = filter.prim_event_index;
  return json;
}

void parse(InputParser<clexmonte::EventFilterGroup> &parser) {
  auto ptr = std::make_unique<clexmonte::EventFilterGroup>();
  clexmonte::EventFilterGroup &filter = *ptr;
  parser.require(filter.unitcell_index, "unitcell_index");
  parser.require(filter.include_by_default, "include_by_default");
  parser.require(filter.prim_event_index, "prim_event_index");
  if (parser.valid()) {
    parser.value = std::move(ptr);
  }
}

void from_json(clexmonte::EventFilterGroup &filter, jsonParser const &json) {
  InputParser<clexmonte::EventFilterGroup> parser{json};
  std::stringstream ss;
  ss << "Error: Invalid clexmonte::EventFilterGroup object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  filter = std::move(*parser.value);
}

}  // namespace CASM
