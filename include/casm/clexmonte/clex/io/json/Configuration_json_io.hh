#ifndef CASM_clexmonte_clex_Configuration_json_io
#define CASM_clexmonte_clex_Configuration_json_io

namespace CASM {

class jsonParser;
template <typename T>
T from_json(jsonParser const &);

namespace clexmonte {
struct Configuration;
}

/// \brief Write clexmonte::Configuration to JSON
jsonParser &to_json(clexmonte::Configuration const &configuration,
                    jsonParser &json);

/// \brief Read clexmonte::Configuration from JSON
template <>
clexmonte::Configuration from_json<clexmonte::Configuration>(
    jsonParser const &json);

/// \brief Read clexmonte::Configuration from JSON
void from_json(clexmonte::Configuration &configuration, jsonParser const &json);

}  // namespace CASM

#endif
