#ifndef CASM_clexmonte_misc_parse_array
#define CASM_clexmonte_misc_parse_array

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

template <typename T, typename... Args>
void parse_array(InputParser<std::vector<T>> &parser, Args &&...args) {
  parser.value = std::make_unique<std::vector<T>>();
  jsonParser const &json = parser.self;
  if (!json.is_array()) {
    parser.error.insert("Error: Expected a JSON array");
    return;
  }

  // for each array element
  Index i = 0;
  for (auto it = json.begin(); it != json.end(); ++it) {
    fs::path relpath = std::to_string(i);
    auto subparser =
        parser.template subparse<T>(relpath, std::forward<Args>(args)...);

    if (subparser->valid()) {
      parser.value->emplace_back(std::move(*(subparser->value)));
    } else {
      return;
    }
    ++i;
  }
}

}  // namespace CASM

#endif
