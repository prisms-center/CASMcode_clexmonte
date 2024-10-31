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
  StateData(std::shared_ptr<system_type> _system, state_type const *_state,
            monte::OccLocation const *_occ_location);

  /// System data (not null)
  std::shared_ptr<system_type> system;

  /// Current state (not null)
  state_type const *state;

  /// Occupant tracker (may be null)
  ///
  /// - This should be the OccLocation being used by the MonteCalculator
  /// - Allows pointing to an OccLocation owned elsewhere
  monte::OccLocation const *occ_location;

  /// Occupant tracker (may be null)
  ///
  /// - Allows constructing a new occupant tracker
  std::shared_ptr<monte::OccLocation> owned_occ_location;

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
  std::map<std::string,
           std::pair<std::shared_ptr<clexulator::MultiClusterExpansion>,
                     std::map<std::string, Index>>>
      multiclex;

  /// Local cluster expansion calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::LocalClusterExpansion>>
      local_clex;

  /// Multi- Local cluster expansion calculators, set for current state
  std::map<std::string,
           std::pair<std::shared_ptr<clexulator::MultiLocalClusterExpansion>,
                     std::map<std::string, Index>>>
      local_multiclex;

  /// Order parameter calculators, set for current state
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;

  /// Local orbit composition calculators, set for current state
  std::map<std::string, std::shared_ptr<LocalOrbitCompositionCalculator>>
      local_orbit_composition_calculators;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
