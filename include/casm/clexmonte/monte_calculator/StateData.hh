#ifndef CASM_clexmonte_StateData
#define CASM_clexmonte_StateData

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/RandomNumberGenerator.hh"

namespace CASM {
namespace clexmonte {

struct StateData {
  StateData(std::shared_ptr<system_type> _system, bool _update_species,
            state_type const *_state, monte::OccLocation const *_occ_location);

  /// System data (not null)
  std::shared_ptr<system_type> system;

  /// Update species in monte::OccLocation tracker?
  bool update_species;

  /// Current state (not null)
  state_type const *state;

  /// Occupant tracker (not null)
  monte::OccLocation const *occ_location;

  /// Current supercell, depends on current state
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Number of unit cells, depends on current state
  Index n_unitcells;

  /// Index conversions, depends on current state (not null)
  monte::Conversions const *convert;

  /// CASM::monte correlation calculators - calculate all correlations
  std::map<std::string, std::shared_ptr<clexulator::Correlations>> corr;

  /// CASM::monte local correlation calculators - calculate all correlations
  std::map<std::string, std::shared_ptr<clexulator::LocalCorrelations>>
      local_corr;

  /// Cluster expansion calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::ClusterExpansion>> clex;

  /// Multi- Cluster expansion calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::MultiClusterExpansion>>
      multiclex;

  /// Local cluster expansion calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::LocalClusterExpansion>>
      local_clex;

  /// Multi- Local cluster expansion calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::MultiLocalClusterExpansion>>
      local_multiclex;

  /// Order parameter calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
