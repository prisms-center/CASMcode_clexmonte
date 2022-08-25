#ifndef CASM_clexmonte_system_OccSystem
#define CASM_clexmonte_system_OccSystem

#include "casm/clexmonte/clex/ClexData.hh"
#include "casm/clexmonte/misc/Matrix3lCompare.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/DoFSpace.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/OrderParameter.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"

namespace CASM {

namespace monte {
template <typename _ConfigType>
struct State;
}

namespace clexmonte {

struct Configuration;
struct OccSystem;
struct OccSystemSupercellData;

/// \brief Info on local cluster expansion basis sets
struct EquivalentsInfo {
  std::vector<clust::IntegralCluster> phenomenal_clusters;
  std::vector<Index> equivalent_generating_op_indices;
  std::vector<xtal::SymOp> equivalent_generating_ops;
};

/// \brief KMC event data
struct OccEventTypeDescription {
  std::vector<occ_events::OccEvent> events;
  std::string kra_clex_name;
  std::string freq_clex_name;
};

/// \brief Data structure for holding Monte Carlo calculation data and methods
///     that should only exist once, and should be accessible by
///     sampling functions - occupation DoF
///
/// Notes:
/// - Use the standalone `get_supercell_data` helper methods to get
///   supercell-specific canonical Monte Carlo calculation data, constructing
///   it as necessary
/// - Use the standalone `get_formation_energy_clex` helper method to get
///   supercell-specific clexulator::ClusterExpansion instance for a given
///   state, constructing it as necessary
struct OccSystem {
  /// \brief Constructor
  OccSystem(std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
            composition::CompositionConverter const &_composition_converter);

  /// Primitive crystal structure and allowed degrees of freedom (DoF)
  std::shared_ptr<config::Prim const> prim;

  // --- TODO: remove duplicates

  /// Primitive crystal structure and allowed degrees of freedom (DoF)
  std::shared_ptr<xtal::BasicStructure const> shared_prim;

  /// Prim global DoF info -- used for DoF value basis conversions
  std::map<DoFKey, xtal::DoFSet> global_dof_info;

  /// Prim local DoF info -- used for DoF value basis conversions
  std::map<DoFKey, std::vector<xtal::SiteDoFSet>> local_dof_info;

  // ---

  /// Composition axes and parametric composition conversions functor
  composition::CompositionConverter composition_converter;

  /// Composition calculation functor
  composition::CompositionCalculator composition_calculator;

  /// Prim neighbor list
  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;

  /// Cluster expansion basis sets
  ///
  /// Maps basis set name -> Clexulator
  std::map<std::string, std::shared_ptr<clexulator::Clexulator>> basis_sets;

  /// Formation energy calculation data and methods. Contains:
  /// - clexulator::PrimNeighborList,
  /// - clexulator::Clexulator,
  /// - clexulator::SparseCoefficients
  std::map<std::string, ClexData> clex_data;

  /// Local cluster expansion basis sets
  ///
  /// Maps local basis set name -> Clexulator
  std::map<std::string, std::shared_ptr<std::vector<clexulator::Clexulator>>>
      local_basis_sets;

  /// Local cluster expansion basis set equivalence info
  ///
  /// Maps local basis set name -> EquivalentsInfo
  std::map<std::string, EquivalentsInfo> equivalents_info;

  /// Local cluster expansion data and methods. Contains:
  /// - clexulator::PrimNeighborList,
  /// - clexulator::Clexulator,
  /// - clexulator::SparseCoefficients
  /// - std::vector<xtal::SymOp> equivalent_generating_ops
  std::map<std::string, LocalClexData> local_clex_data;

  /// DoFSpace that define order parameters
  std::map<std::string, clexulator::DoFSpace> order_parameter_definitions;

  /// Supercell specific formation energy calculation data and methods (using
  /// transformation_matrix_to_super as key). Contains:
  /// -  clexulator::SuperNeighborList,
  /// -  clexulator::ClusterExpansion
  std::map<Eigen::Matrix3l, OccSystemSupercellData, Matrix3lCompare>
      supercell_data;

  /// KMC events index definitions
  std::shared_ptr<occ_events::OccSystem> event_system;

  /// KMC events
  std::map<std::string, OccEventTypeData> event_types;

};

/// \brief Data structure for holding supercell-specific Monte Carlo calculation
///     data and methods that should only exist once, and should be accessible
///     by sampling functions - occupation DoF
struct OccSystemSupercellData {
  /// \brief Constructor
  OccSystemSupercellData(OccSystem const &system_data,
                         Eigen::Matrix3l const &transformation_matrix_to_super);

  /// Performs index conversions in supercell
  monte::Conversions convert;

  /// List of unique pairs of (asymmetric unit index, species index)
  monte::OccCandidateList occ_candidate_list;

  /// SuperNeighborList, used for evaluating correlations in a particular
  /// supercell
  std::shared_ptr<clexulator::SuperNeighborList> supercell_neighbor_list;

  /// CASM::monte compatible formation energy calculator. Contains:
  /// -  clexulator::Correlations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::ClusterExpansion>> clex;

  /// CASM::monte compatible formation energy calculator. Contains:
  /// -  clexulator::LocalCorrelations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::LocalClusterExpansion>>
      local_clex;

  /// Order parameter calculators
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;
};

// ---
// The following are used to construct a common interface between "System"
// data, in this case OccSystem, and templated CASM::clexmonte methods such as
// sampling function factory methods
// ---

/// \brief Helper to get std::shared_ptr<xtal::BasicStructure const>
std::shared_ptr<xtal::BasicStructure const> const &get_shared_prim(
    OccSystem &data);

/// \brief Helper to get composition::CompositionConverter
composition::CompositionConverter const &get_composition_converter(
    OccSystem &data);

/// \brief Helper to get composition::CompositionCalculator
composition::CompositionCalculator const &get_composition_calculator(
    OccSystem &data);

/// \brief Helper to make the default configuration in prim basis
Configuration make_default_configuration(
    OccSystem const &data,
    Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    OccSystem const &data,
    Configuration const &configuration_in_standard_basis);

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    OccSystem const &data, Configuration const &configuration_in_prim_basis);

/// \brief Helper to get the ClexData for formation energy
ClexData &get_formation_energy_clex_data(OccSystem &data);

/// \brief Helper to get the ClexData for formation energy
ClexData &get_clex_data(OccSystem &data, std::string const &key);

/// \brief Helper to get LocalClexData
LocalClexData &get_local_clex_data(OccSystem &data, std::string const &key);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular configuration, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_formation_energy_clex(
    OccSystem &data, Configuration const &configuration);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_formation_energy_clex(
    OccSystem &data, monte::State<Configuration> const &state);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular configuration, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_clex(
    OccSystem &data, std::string const &key, Configuration const &configuration);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_clex(
    OccSystem &data, std::string const &key, monte::State<Configuration> const &state);

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular configuration, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    OccSystem &data, std::string const &key,
    Configuration const &configuration);

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    OccSystem &data, std::string const &key,
    monte::State<Configuration> const &state);

/// \brief Helper to get the correct order parameter calculators for a
///     particular configuration, constructing as necessary
std::map<std::string, std::shared_ptr<clexulator::OrderParameter>> const &
get_order_parameters(OccSystem &data, Configuration const &configuration);

/// \brief Helper to get the correct order parameter calculators for a
///     particular state's supercell, constructing as necessary
std::map<std::string, std::shared_ptr<clexulator::OrderParameter>> const &
get_order_parameters(OccSystem &data, monte::State<Configuration> const &state);

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(
    OccSystem &data, monte::State<Configuration> const &state);

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(
    OccSystem &data, monte::State<Configuration> const &state);

}  // namespace clexmonte
}  // namespace CASM

#endif
