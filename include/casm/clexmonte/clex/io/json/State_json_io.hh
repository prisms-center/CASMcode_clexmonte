#ifndef CASM_clexmonte_clex_State_json_io
#define CASM_clexmonte_clex_State_json_io

namespace CASM {

class jsonParser;
template <typename T>
T from_json(jsonParser const &);

namespace clexmonte {
struct Configuration;
}

namespace monte {
template <typename _ConfigType>
struct State;
}

/// \brief Write monte::State<clexmonte::Configuration> to JSON
jsonParser &to_json(monte::State<clexmonte::Configuration> const &state,
                    jsonParser &json);

/// \brief Read monte::State<clexmonte::Configuration> from JSON
template <>
monte::State<clexmonte::Configuration>
from_json<monte::State<clexmonte::Configuration>>(jsonParser const &json);

/// \brief Read monte::State<clexmonte::Configuration> from JSON
void from_json(monte::State<clexmonte::Configuration> &state,
               jsonParser const &json);

}  // namespace CASM

#endif
