#ifndef CASM_clexmonte_misc_to_json
#define CASM_clexmonte_misc_to_json

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

template <typename T>
jsonParser qto_json(T const &t) {
  jsonParser json;
  to_json(t, json);
  return json;
}

}  // namespace CASM

#endif
