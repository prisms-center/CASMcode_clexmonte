#ifndef CASM_clexmonte_system_System
#define CASM_clexmonte_system_System

#include "casm/clexmonte/misc/Matrix3lCompare.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/DoFSpace.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/OrderParameter.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
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
struct System;
struct SupercellSystemData;

/// \brief Data structure for holding Monte Carlo calculation data and methods
///     that should only exist once, and should be accessible by
///     sampling functions - occupation DoF
///
/// Notes:
/// - Use the standalone `get_supercell_data` helper methods to get
///   supercell-specific canonical Monte Carlo calculation data, constructing
///   it as necessary
/// - Use the standalone `get_clex` helper method to get
///   supercell-specific clexulator::ClusterExpansion instance for a given
///   state, constructing it as necessary
struct System {
  /// \brief Constructor
  System(std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
         composition::CompositionConverter const &_composition_converter);

  // --- Crystal structure

  /// Primitive crystal structure and allowed degrees of freedom (DoF)
  std::shared_ptr<config::Prim const> prim;

  // --- Composition

  /// Composition axes and parametric composition conversions functor
  composition::CompositionConverter composition_converter;

  /// Composition calculation functor
  composition::CompositionCalculator composition_calculator;

  // --- Order parameters

  /// DoFSpace that define order parameters
  std::map<std::string, clexulator::DoFSpace> order_parameter_definitions;

  // --- Cluster expansions

  /// Prim neighbor list
  ///
  /// Notes:
  /// - Make sure to use this->prim_neighbor_list to construct
  ///   Clexulators stored in this->basis_sets and this->local_basis_sets
  ///   so that SupercellSystemData can properly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;

  /// Cluster expansion basis sets
  ///
  /// Notes:
  /// - Maps basis set name -> Clexulator
  /// - Make sure to use this->prim_neighbor_list to construct the
  ///   Clexulator so that SupercellSystemData can propertly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::map<std::string, std::shared_ptr<clexulator::Clexulator>> basis_sets;

  /// Data used to construct clexulator::ClusterExpansion. Contains:
  /// - basis_set_name
  /// - clexulator::SparseCoefficients
  std::map<std::string, ClexData> clex_data;

  /// Data used to construct clexulator::MultiClusterExpansion. Contains:
  /// - basis_set_name
  /// - std::vector<clexulator::SparseCoefficients>
  std::map<std::string, MultiClexData> multiclex_data;

  // --- Local cluster expansions

  /// Local cluster expansion basis sets
  ///
  /// Notes:
  /// - Maps local basis set name -> std::vector<clexulator::Clexulator>
  /// - Make sure to use this->prim_neighbor_list to construct the
  ///   Clexulator so that SupercellSystemData can propertly construct
  ///   the SuperNeighborList needed to evaluate correlations
  std::map<std::string, std::shared_ptr<std::vector<clexulator::Clexulator>>>
      local_basis_sets;

  /// Local cluster expansion basis set equivalence info
  ///
  /// Maps local basis set name -> EquivalentsInfo
  std::map<std::string, EquivalentsInfo> equivalents_info;

  /// Data used to construct clexulator::LocalClusterExpansion. Contains:
  /// - local_basis_set_name
  /// - clexulator::SparseCoefficients
  std::map<std::string, LocalClexData> local_clex_data;

  /// Data used to construct clexulator::LocalMultiClusterExpansion. Contains:
  /// - local_basis_set_name
  /// - std::vector<clexulator::SparseCoefficients>
  std::map<std::string, LocalMultiClexData> local_multiclex_data;

  // --- KMC events

  /// KMC events index definitions
  std::shared_ptr<occ_events::OccSystem> event_system;

  /// KMC event symgroup representation
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;

  /// KMC events
  std::map<std::string, OccEventTypeData> event_type_data;

  // --- Supercells

  /// Supercell specific formation energy calculation data and methods (using
  /// transformation_matrix_to_super as key).
  std::map<Eigen::Matrix3l, SupercellSystemData, Matrix3lCompare>
      supercell_data;
};

/// \brief Data structure for holding supercell-specific Monte Carlo calculation
///     data and methods that should only exist once, and should be accessible
///     by sampling functions - occupation DoF
struct SupercellSystemData {
  /// \brief Constructor
  SupercellSystemData(System const &system,
                      Eigen::Matrix3l const &transformation_matrix_to_super);

  /// Number of unit cells in the supercell
  Index n_unitcells;

  // --- Index conversions and occupation tracking

  /// Performs index conversions in supercell
  monte::Conversions convert;

  /// List of unique pairs of (asymmetric unit index, species index)
  monte::OccCandidateList occ_candidate_list;

  // --- Cluster expansion

  /// SuperNeighborList, used for evaluating correlations in a particular
  /// supercell
  std::shared_ptr<clexulator::SuperNeighborList> supercell_neighbor_list;

  // --- Order parameter

  /// Order parameter calculators
  std::map<std::string, std::shared_ptr<clexulator::OrderParameter>>
      order_parameters;

  // --- Cluster expansion

  /// CASM::monte compatible cluster expansion calculators. Contains:
  /// -  clexulator::Correlations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::ClusterExpansion>> clex;

  /// CASM::monte compatible cluster expansion calculators. Contains:
  /// -  clexulator::Correlations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::MultiClusterExpansion>>
      multiclex;

  /// CASM::monte compatible local cluster expansion calculators. Contains:
  /// -  clexulator::LocalCorrelations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::LocalClusterExpansion>>
      local_clex;

  /// CASM::monte compatible local cluster expansion calculators. Contains:
  /// -  clexulator::LocalCorrelations
  /// -  clexulator::SparseCoefficients
  std::map<std::string, std::shared_ptr<clexulator::MultiLocalClusterExpansion>>
      local_multiclex;
};

// ---
// The following are used to construct a common interface between "System"
// data, in this case System, and templated CASM::clexmonte methods such as
// sampling function factory methods
// ---

/// \brief Helper to get std::shared_ptr<config::Prim const>
std::shared_ptr<config::Prim const> const &get_prim_info(System const &system);

/// \brief Helper to get std::shared_ptr<xtal::BasicStructure const>
std::shared_ptr<xtal::BasicStructure const> const &get_prim_basicstructure(
    System const &system);

/// \brief Helper to get prim basis
std::vector<xtal::Site> const &get_basis(System const &system);

/// \brief Helper to get basis size
Index get_basis_size(System const &system);

/// \brief Helper to get composition::CompositionConverter
composition::CompositionConverter const &get_composition_converter(
    System const &system);

/// \brief Helper to get composition::CompositionCalculator
composition::CompositionCalculator const &get_composition_calculator(
    System const &system);

/// \brief Helper to make the default configuration in prim basis
Configuration make_default_configuration(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system, Configuration const &configuration_in_standard_basis);

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis);

/// \brief Helper to make the default configuration in prim basis
monte::State<Configuration> make_default_state(
    System const &system,
    Eigen::Matrix3l const &transformation_matrix_to_super);

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    System const &system, Configuration const &configuration_in_standard_basis);

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    System const &system, Configuration const &configuration_in_prim_basis);

/// \brief Helper to get the Clexulator
std::shared_ptr<clexulator::Clexulator> get_basis_set(System const &system,
                                                      std::string const &key);

/// \brief Helper to get the local Clexulator
std::shared_ptr<std::vector<clexulator::Clexulator>> get_local_basis_set(
    System const &system, std::string const &key);

/// \brief Helper to get ClexData
ClexData const &get_clex_data(System const &system, std::string const &key);

/// \brief Helper to get MultiClexData
MultiClexData const &get_multiclex_data(System const &system,
                                        std::string const &key);

/// \brief Helper to get LocalClexData
LocalClexData const &get_local_clex_data(System const &system,
                                         std::string const &key);

/// \brief Helper to get LocalMultiClexData
LocalMultiClexData const &get_local_multiclex_data(System const &system,
                                                   std::string const &key);

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalClexData const &local_clex_data,
    Index equivalent_index);

/// \brief Construct impact tables
std::set<xtal::UnitCellCoord> get_required_update_neighborhood(
    System const &system, LocalMultiClexData const &local_multiclex_data,
    Index equivalent_index);

/// \brief KMC events index definitions
std::shared_ptr<occ_events::OccSystem> get_event_system(System const &system);

/// \brief KMC events
std::map<std::string, OccEventTypeData> const &get_event_type_data(
    System const &system);

/// \brief KMC events
OccEventTypeData const &get_event_type_data(System const &system,
                                            std::string const &key);

// --- Supercell-specific

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_clex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key);

/// \brief Helper to get the correct clexulator::MultiClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::MultiClusterExpansion> get_multiclex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key);

/// \brief Helper to get the correct clexulator::LocalClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::LocalClusterExpansion> get_local_clex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key);

/// \brief Helper to get the correct clexulator::MultiLocalClusterExpansion for
/// a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::MultiLocalClusterExpansion> get_local_multiclex(
    System &system, monte::State<Configuration> const &state,
    std::string const &key);

/// \brief Helper to get the correct order parameter calculators for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::OrderParameter> get_order_parameter(
    System &system, monte::State<Configuration> const &state,
    std::string const &key);

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(
    System &system, monte::State<Configuration> const &state);

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(
    System &system, monte::State<Configuration> const &state);

}  // namespace clexmonte
}  // namespace CASM

#endif