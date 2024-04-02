#ifndef CASM_clexmonte_RunData
#define CASM_clexmonte_RunData

#include <optional>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/ValueMap.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

struct RunData {
  std::optional<state_type> initial_state;
  std::optional<state_type> final_state;
  monte::ValueMap conditions;
  Eigen::Matrix3l transformation_matrix_to_super;
  Index n_unitcells;
};

struct RunDataOutputParams {
  /// \brief Save all initial_state in completed_runs
  bool do_save_all_initial_states = false;

  /// \brief Save all final_state in completed_runs
  bool do_save_all_final_states = false;

  /// \brief Save last final_state in completed_runs
  bool do_save_last_final_state = true;

  /// \brief Write saved initial_state to completed_runs.json
  bool write_initial_states = false;

  /// \brief Write saved final_state to completed_runs.json
  bool write_final_states = false;

  /// \brief Location to save completed_runs.json if not empty
  fs::path output_dir;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
