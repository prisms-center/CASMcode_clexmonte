#ifndef CASM_clexmonte_run_ResultsIO_json_io
#define CASM_clexmonte_run_ResultsIO_json_io

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/polymorphic_method_json_io.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Construct ResultsIO from JSON
void parse(InputParser<results_io_type> &parser,
           MethodParserMap<results_io_type> const &results_io_methods);

/// \brief Construct jsonResultsIO from JSON
void parse(InputParser<monte::jsonResultsIO<config_type>> &parser,
           std::map<std::string, state_sampling_function_type> const
               &sampling_functions,
           std::map<std::string, results_analysis_function_type> const
               &analysis_functions);

}  // namespace clexmonte
}  // namespace CASM

#endif
