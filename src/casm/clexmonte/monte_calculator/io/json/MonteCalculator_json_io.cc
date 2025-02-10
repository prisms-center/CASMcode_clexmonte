#include "casm/clexmonte/monte_calculator/io/json/MonteCalculator_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {
namespace clexmonte {

namespace {

fs::path resolve_path(fs::path p, std::vector<fs::path> search_path) {
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

}  // namespace

/// \brief Parse MonteCalculator from JSON
///
/// \param parser An InputParser, as genrated by
///     `InputParser::subparse<std::shared_ptr<MonteCalculator>>` or one of the
///     other `subparse` methods.
/// \param system System data
/// \param params Calculation method parameters, as specified by the
///     particular calculation type
/// \param engine A random number engine. If null, a new one seeded by
///     std::random_device is constructed.
/// \param search_path Paths besides the current working directory to look for
///     a MonteCalculator source file.
///
/// Expected JSON format:
///   source: string (required)
///     Path to a MonteCalculator source file implementing a Monte Carlo
///     calculator.
///
///   compile_options: (optional)
///     Options used to compile the MonteCalculator source file, if it is not
///     yet compiled. Example:
///     "g++ -O3 -Wall -fPIC --std=c++17 -I/path/to/include"
///
///   so_options: (optional)
///     Options used to compile the MonteCalculator shared object file, if it
///     is not yet compiled. Example:
///     "g++ -shared -L/path/to/lib -lcasm_clexmonte "
void parse(InputParser<std::shared_ptr<MonteCalculator>> &parser,
           std::shared_ptr<System> &system, jsonParser const &params,
           std::shared_ptr<MonteCalculator::engine_type> engine,
           std::vector<fs::path> search_path) {
  // parse "source"
  std::string _calculator_src;
  parser.require(_calculator_src, "source");
  fs::path calculator_src = resolve_path(_calculator_src, search_path);
  if (!fs::exists(calculator_src)) {
    parser.insert_error("source", "Error: \"source\" file does not exist.");
  }

  if (!parser.valid()) {
    return;
  }

  // - name of MonteCalculator source file (excluding .cc extension)
  std::string calculator_name = calculator_src.stem();

  // - directory where the MonteCalculator source file is found
  fs::path calculator_dirpath = calculator_src.parent_path();

  // - set MonteCalculator compilation options
  //   ex: g++ -O3 -Wall -fPIC --std=c++17 -I/path/to/include
  std::string default_calculator_compile_options =
      //
      // uses $CASM_CXX, else default="g++"
      RuntimeLibrary::default_cxx().first + " " +
      //
      // uses $CASM_CXXFLAGS, else default="-O3 -Wall -fPIC --std=c++17"
      RuntimeLibrary::default_cxxflags().first + " " +
      //
      // uses -I$CASM_INCLUDEDIR,
      //   else -I$CASM_PREFIX/include,
      //   else tries to find "ccasm" or "casm" executable on PATH and looks
      //     for standard include paths relative from there,
      //   else fails with "/not/found"
      include_path(RuntimeLibrary::default_casm_includedir().first);

  std::string calculator_compile_options;
  parser.optional_else(calculator_compile_options, "compile_options",
                       default_calculator_compile_options);

  // - set MonteCalculator shared object compilation options
  //   ex: g++ -shared -L/path/to/lib -lcasm_global -lcasm_crystallography
  //     -lcasm_clexulator -lcasm_monte -lcasm_clexmonte
  std::string default_calculator_so_options =
      //
      // uses $CASM_CXX, else default="g++"
      RuntimeLibrary::default_cxx().first + " " +
      //
      // uses $CASM_SOFLAGS, else default="-shared"
      RuntimeLibrary::default_soflags().first + " " +
      //
      // uses -L$CASM_LIBDIR,
      //   else -L$CASM_PREFIX/lib,
      //   else tries to find "ccasm" or "casm" executables on PATH and looks
      //     for libcasm at standard relative paths from there,
      //   else fails with "-L/not/found"
      link_path(RuntimeLibrary::default_casm_libdir().first) + " " +
      //
      // requires libcasm_clexmonte:
      "-lcasm_clexmonte ";

  std::string calculator_so_options;
  parser.optional_else(calculator_so_options, "so_options",
                       default_calculator_so_options);

  if (parser.valid()) {
    parser.value = std::make_unique<std::shared_ptr<MonteCalculator>>(
        make_monte_calculator_from_source(
            calculator_dirpath, calculator_name, system, params, engine,
            calculator_compile_options, calculator_so_options));
  }
}

}  // namespace clexmonte
}  // namespace CASM
