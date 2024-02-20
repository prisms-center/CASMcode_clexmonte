#ifndef CASM_clexmonte_StateGenerator
#define CASM_clexmonte_StateGenerator

#include <vector>

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/run/io/json/RunData_json_io.hh"
#include "casm/monte/definitions.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

/// \brief A StateGenerator generates initial states for a series of Monte
///     Carlo calculations
///
/// Notes:
/// - The template parameter _RunInfoType is specified by a particular Monte
///   Carlo method implementation.
/// - _RunInfoType allows customization of what
///   information is provided to a particular state generation method. In the
///   basic case, it will be the final state for each run. Templating allows
///   support for more complex cases where the next state could be generated
///   based on the sampled data collected during previous runs.
class StateGenerator {
 public:
  virtual ~StateGenerator() {}

  /// \brief Check if calculations are complete, using info from all finished
  ///     runs
  virtual bool is_complete() = 0;

  /// \brief Generate the next initial state, using info from all finished
  ///     runs
  virtual state_type next_state() = 0;

  /// \brief Push back data for a completed run
  virtual void push_back(RunData const &run_data) = 0;

  virtual Index n_completed_runs() const = 0;

  virtual std::vector<RunData> const &completed_runs() const = 0;

  virtual void read_completed_runs() = 0;

  virtual void write_completed_runs() const = 0;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
