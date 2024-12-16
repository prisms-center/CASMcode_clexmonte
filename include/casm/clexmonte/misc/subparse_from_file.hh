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
  std::string filepath;
  parser.require(filepath, option);

  fs::path resolved_path = resolve_path(filepath, search_path);

  if (!fs::exists(resolved_path)) {
    parser.insert_error(option, "Error: file not found.");
    jsonParser json;
    auto subparser = std::make_shared<InputParser<RequiredType>>(
        json, std::forward<Args>(args)...);
    return subparser;
  }
  jsonParser json{resolved_path};
  auto subparser = std::make_shared<InputParser<RequiredType>>(
      json, std::forward<Args>(args)...);

  auto &log = CASM::log();
  if (!subparser->valid()) {
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
      log.indent() << report << std::endl << std::endl;
    }

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
  }
  if (subparser->all_warnings().size()) {
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

    parser.insert_warning(
        option, "Warning: warnings for file: " + resolved_path.string());
    for (auto const &warning : subparser->all_warnings()) {
      for (auto const &msg : warning.second) {
        parser.insert_warning(option,
                              "@(/" + warning.first.string() + "): " + msg);
      }
    }
  }

  subparser->type_name = CASM::type_name<RequiredType>();
  return subparser;
}

}  // namespace CASM

#endif
