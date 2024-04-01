#ifndef CASM_clexmonte_clex_State_json_io
#define CASM_clexmonte_clex_State_json_io

#include "casm/clexmonte/definitions.hh"

namespace CASM {

namespace config {
class SupercellSet;
}

class jsonParser;
template <typename T>
class InputParser;
template <typename T>
T from_json(jsonParser const &);
template <typename T>
struct jsonConstructor;

/// \brief Write monte::State<clexmonte::Configuration> to JSON
jsonParser &to_json(monte::State<clexmonte::Configuration> const &state,
                    jsonParser &json, bool write_prim_basis = false);

void parse(InputParser<monte::State<clexmonte::Configuration>> &parser,
           config::SupercellSet &supercells);

/// \brief Read monte::State<clexmonte::Configuration> from JSON
template <>
struct jsonConstructor<monte::State<clexmonte::Configuration>> {
  static monte::State<clexmonte::Configuration> from_json(
      jsonParser const &json, config::SupercellSet &supercells);
};

/// \brief Read monte::State<clexmonte::Configuration> from JSON
void from_json(monte::State<clexmonte::Configuration> &state,
               jsonParser const &json, config::SupercellSet &supercells);

}  // namespace CASM

#endif
