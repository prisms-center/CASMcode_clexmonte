#ifndef CASM_clexmonte_ConfigGenerator
#define CASM_clexmonte_ConfigGenerator

#include <vector>

#include "casm/clexmonte/definitions.hh"

namespace CASM {
namespace clexmonte {

/// \brief A ConfigGenerator generates a configuration given a set of
///     conditions and results from previous runs
///
/// Notes:
/// - The template parameter _RunInfoType is specified by a particular Monte
///   Carlo method implementation.
/// - _RunInfoType allows customization of what
///   information is provided to a particular configuration generation method.
///   In the basic case, it will be the final state for each run. Templating
///   allows support for more complex cases where the next state could be
///   generated based on the sampled data collected during previous runs.
class ConfigGenerator {
 public:
  virtual ~ConfigGenerator() {}

  /// \brief Generate a configuration, using information from a set of
  /// conditions and info from previous runs
  virtual config_type operator()(
      monte::ValueMap const &conditions,
      std::vector<RunData> const &completed_runs) = 0;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
