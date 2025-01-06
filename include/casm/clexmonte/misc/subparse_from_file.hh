#ifndef CASM_clexmonte_misc_subparse_from_file
#define CASM_clexmonte_misc_subparse_from_file

#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

inline fs::path resolve_path(fs::path p, std::vector<fs::path> search_path) {
  fs::path resolved_path;
  if (fs::exists(p)) {
    return p;
  } else {
    for (fs::path root : search_path) {
      if (fs::exists(root / p)) {
        return root / p;
      }
    }
  }
  return p;
}

/// Run an InputParser on the JSON file with path given by the option,
///     collecting errors and warnings
///
/// \param parser The InputParser
/// \param option The option that gives a file path
/// \param search_path A vector of paths to use as the root to resolve
///     the file path given by `option`, if that file path is a relative
///     path.
/// \param args Additional args to pass to subparser
///
template <typename RequiredType, typename T, typename... Args>
std::shared_ptr<InputParser<RequiredType>> subparse_from_file(
    InputParser<T> &parser, fs::path option,
    std::vector<fs::path> search_path = {}, Args &&...args) {
  // Logging for debug verbosity only
  // (except for errors and warnings noted below)
  auto &log = CASM::log();
  log.begin_section<Log::verbose>();
  log.increase_indent();
  log.indent() << "- Attempting to subparse from file..." << std::endl;

  jsonParser null_json;
  auto null_subparser = std::make_shared<InputParser<RequiredType>>(
      null_json, std::forward<Args>(args)...);

  auto it = parser.self.find_at(option);
  if (it == parser.self.end()) {
    std::stringstream msg;
    msg << "Error: missing required option '" << option.string() << "'.";
    parser.insert_error(option, msg.str());

    return null_subparser;
  }

  std::string filepath;
  parser.require(filepath, option);
  log.indent() << "- filepath: " << filepath << std::endl;

  fs::path resolved_path = resolve_path(filepath, search_path);
  log.indent() << "- resolved_path: " << resolved_path << std::endl;

  if (!fs::exists(resolved_path)) {
    log.indent() << "- resolved_path not found" << std::endl;
    log.decrease_indent();
    log.end_section();

    parser.insert_error(option, "Error: file not found.");
    return null_subparser;
  }

  log.indent() << "- Reading file..." << std::endl;
  jsonParser json{resolved_path};
  log.indent() << "- Subparsing..." << std::endl;
  auto subparser = std::make_shared<InputParser<RequiredType>>(
      json, std::forward<Args>(args)...);

  if (!subparser->valid()) {
    log.indent() << "- Subparsing: failed" << std::endl << std::endl;

    // Always log this:
    log.begin_section<Log::none>();
    log << std::endl;
    log << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    log << "~~~  Failed parsing input file  ~~~" << std::endl;
    log << std::endl;
    log << "file: " << resolved_path << std::endl;
    log << std::endl;
    print_errors(*subparser, log, "Error Summary");
    log << std::endl;
    if (subparser->all_warnings().size()) {
      print_warnings(*subparser, log, "Warning Summary");
      log << std::endl;
    }
    if (json.is_obj()) {
      jsonParser report = make_report(*subparser);
      log << report << std::endl << std::endl;
    }
    log.end_section();

    parser.insert_error(
        option, "Error: Failed to parse file: " + resolved_path.string());
    for (auto const &error : subparser->all_errors()) {
      for (auto const &msg : error.second) {
        parser.insert_error(option, "@(/" + error.first.string() + "): " + msg);
      }
    }
    for (auto const &warning : subparser->all_warnings()) {
      for (auto const &msg : warning.second) {
        parser.insert_warning(option,
                              "@(/" + warning.first.string() + "): " + msg);
      }
    }
  } else {
    log.indent() << "- Subparsing: succeeded" << std::endl;
  }

  if (subparser->all_warnings().size()) {
    // Always log this:
    log.begin_section<Log::none>();
    log << std::endl;
    log << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    log << "~~~ Warnings parsing input file ~~~" << std::endl;
    log << std::endl;
    log << "file: " << resolved_path << std::endl;
    log << std::endl;
    print_warnings(*subparser, log, "Warning Summary");
    log << std::endl;
    if (json.is_obj()) {
      jsonParser report = make_report(*subparser);
      log.indent() << report << std::endl << std::endl;
    }
    log.end_section();

    parser.insert_warning(
        option, "Warning: warnings for file: " + resolved_path.string());
    for (auto const &warning : subparser->all_warnings()) {
      for (auto const &msg : warning.second) {
        parser.insert_warning(option,
                              "@(/" + warning.first.string() + "): " + msg);
      }
    }
  }

  log.indent() << "- Subparsing finished..." << std::endl;
  log.end_section();
  log.decrease_indent();

  subparser->type_name = CASM::type_name<RequiredType>();
  return subparser;
}

}  // namespace CASM

#endif
