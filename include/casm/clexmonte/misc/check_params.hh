#ifndef CASM_clexmonte_misc_check_params
#define CASM_clexmonte_misc_check_params

#include <optional>
#include <set>
#include <sstream>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
namespace clexmonte {

/// \brief Check if JSON parameters have all required parameters, and otherwise
/// only optional parameters or "comment" parameters which begin with '_`.
///
/// Raises if any required parameters are missing.
/// Prints a warning for any unknown parameters to the CASM::err_log().
///
/// \param params JSON object with parameters
/// \param required_params Required parameters
/// \param optional_params Optional parameters
/// \param location Path of the JSON object being checked, for error messages
inline void check_params(jsonParser const &params,
                         std::set<std::string> const &required_params,
                         std::set<std::string> const &optional_params,
                         std::optional<fs::path> location = std::nullopt) {
  for (auto key : required_params) {
    if (!params.contains(key)) {
      std::stringstream msg;
      fs::path key_path(key);
      if (location) {
        key_path = *location / key_path;
      }
      msg << "Error: Missing required parameter '" << key_path.string() << "'.";
      throw std::runtime_error(msg.str());
    }
  }

  auto &log = CASM::log();
  std::vector<std::string> unknown_parameters;
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.name();
    if (key.empty()) {
      std::stringstream msg;
      msg << "Error: Empty parameter key.";
      throw std::runtime_error(msg.str());
    }
    if (key[0] == '_') {
      continue;
    }
    if (!required_params.count(key) && !optional_params.count(key)) {
      unknown_parameters.push_back(key);
      log.indent() << "Warning: Unknown parameter '" << key << "'."
                   << std::endl;
    }
  }

  if (unknown_parameters.size()) {
    Log &log = CASM::err_log();
    log << std::endl;
    log << "## WARNING: UNKNOWN PARAMETERS ######################\n"
           "#                                                   #\n"
           "# Parameters:                                       #\n";
    for (auto const &key : unknown_parameters) {
      fs::path key_path(key);
      if (location) {
        key_path = *location / key_path;
      }
      log << "  - '" << key_path.string() << "'\n";
    }
    log << "#                                                   #\n"
           "#####################################################\n"
        << std::endl;
  }
}

}  // namespace clexmonte
}  // namespace CASM
#endif
