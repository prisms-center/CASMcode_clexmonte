#ifndef CASM_clexmonte_system_OccSystem
#define CASM_clexmonte_system_OccSystem

#include "casm/clexmonte/clex/ClexData.hh"
#include "casm/clexmonte/misc/Matrix3lCompare.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/NeighborList.hh"
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
            composition::CompositionConverter const &_composition_converter,
            ClexData const &_formation_energy_clex_data);

  /// Primitive crystal structure and allowed degrees of freedom (DoF)
  std::shared_ptr<xtal::BasicStructure const> shared_prim;

  /// Prim global DoF info -- used for DoF value basis conversions
  std::map<DoFKey, xtal::DoFSet> global_dof_info;

  /// Prim local DoF info -- used for DoF value basis conversions
  std::map<DoFKey, std::vector<xtal::SiteDoFSet>> local_dof_info;

  /// Composition axes and parametric composition conversions functor
  composition::CompositionConverter composition_converter;

  /// Composition calculation functor
  composition::CompositionCalculator composition_calculator;

  /// Formation energy calculation data and methods. Contains:
  /// - clexulator::PrimNeighborList,
  /// - clexulator::Clexulator,
  /// - clexulator::SparseCoefficients
  ClexData formation_energy_clex_data;

  /// Supercell specific formation energy calculation data and methods (using
  /// transformation_matrix_to_super as key). Contains:
  /// -  clexulator::SuperNeighborList,
  /// -  clexulator::ClusterExpansion
  std::map<Eigen::Matrix3l, OccSystemSupercellData, Matrix3lCompare>
      supercell_data;
};

/// \brief Data structure for holding supercell-specific Monte Carlo calculation
///     data and methods that should only exist once, and should be accessible
///     by sampling functions - occupation DoF
struct OccSystemSupercellData {
  /// \brief Constructor
  OccSystemSupercellData(OccSystem const &system_data,
                         Eigen::Matrix3l const &transformation_matrix_to_super);

  /// SuperNeighborList, used for evaluating correlations in a particular
  /// supercell
  std::shared_ptr<clexulator::SuperNeighborList> supercell_neighbor_list;

  /// CASM::monte compatible formation energy calculator. Contains:
  /// -  clexulator::Correlations
  /// -  clexulator::SparseCoefficients
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;

  /// Performs index conversions in supercell
  monte::Conversions convert;

  /// List of unique pairs of (asymmetric unit index, species index)
  monte::OccCandidateList occ_candidate_list;
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

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular configuration, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_formation_energy_clex(
    OccSystem &data, Configuration const &configuration);

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state's supercell, constructing as necessary
std::shared_ptr<clexulator::ClusterExpansion> get_formation_energy_clex(
    OccSystem &data, monte::State<Configuration> const &state);

/// \brief Helper to get supercell index conversions
monte::Conversions const &get_index_conversions(
    OccSystem &data, monte::State<Configuration> const &state);

/// \brief Helper to get unique pairs of (asymmetric unit index, species index)
monte::OccCandidateList const &get_occ_candidate_list(
    OccSystem &data, monte::State<Configuration> const &state);

}  // namespace clexmonte
}  // namespace CASM

#endif
