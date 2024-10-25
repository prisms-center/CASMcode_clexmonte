#ifndef CASM_clexmonte_LocalOrbitCompositionCalculator
#define CASM_clexmonte_LocalOrbitCompositionCalculator

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/monte/misc/LexicographicalCompare.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
}

namespace composition {
class CompositionCalculator;
}

namespace clexmonte {

struct PossibleLocalOrbitCompositions {
  /// \brief The maximum number of possible unique local orbit compositions
  Index max_possible_compositions;

  /// \brief Whether all possible unique local orbit compositions have been
  /// found
  bool found_all;

  /// \brief Number of components
  Index n_components;

  /// \brief The number of local orbits composition is calculated for
  Index n_orbits;

  /// \brief  Holds unique num_each_component_by_orbit as column-major
  ///     Eigen::VectorXl, up to max_possible_compositions
  std::set<Eigen::VectorXl, monte::LexicographicalCompare> possible;

  /// \brief Calculate the number of possible unique local orbit compositions
  PossibleLocalOrbitCompositions(
      composition::CompositionCalculator const &comp_calculator,
      std::vector<std::vector<std::set<clust::IntegralCluster>>> const &orbits,
      Index _max_possible_occupations);
};

/// \brief Make a component name, which describes a composition matrix using
/// orbit indices and the component names
std::string make_component_name(
    Eigen::MatrixXi const &num_each_component_by_orbit,
    std::set<int> orbits_to_calculate,
    std::vector<std::string> const &components, bool combine_orbits);

class LocalOrbitCompositionCalculator {
 public:
  //  LocalOrbitCompositionCalculator(std::shared_ptr<system_type> _system,
  //                                  std::string _event_type_name,
  //                                  std::set<int> _orbits_to_calculate);
  //
  //  /// \brief Reset pointer to state currently being calculated
  //  void set(state_type const *state);

  /// \brief Constructor - for a single supercell
  LocalOrbitCompositionCalculator(
      std::vector<std::vector<std::set<clust::IntegralCluster>>> const &_orbits,
      std::set<int> _orbits_to_calculate, bool _combine_orbits,
      std::shared_ptr<clexulator::PrimNeighborList> _prim_nlist,
      std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
      xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
      composition::CompositionCalculator const &_composition_calculator,
      clexulator::ConfigDoFValues const *_dof_values = nullptr);

  /// \brief Reset pointer to configuration currently being calculated
  void set(clexulator::ConfigDoFValues const *dof_values);

  /// \brief Calculate the composition by orbit around an event
  Eigen::MatrixXi const &value(Index unitcell_index, Index equivalent_index);

  /// \brief The orbits to calculate
  std::set<int> const &orbits_to_calculate() const {
    return m_orbits_to_calculate;
  }

  /// \brief Whether orbits are combined
  bool combine_orbits() const { return m_combine_orbits; }

  /// \brief The local orbit neighbor and sublattice indices
  std::vector<std::vector<std::set<std::pair<int, int>>>> const &
  local_orbits_neighbor_indices() const {
    return m_local_orbits_neighbor_indices;
  }

  /// \brief The local orbit sites
  std::vector<std::vector<std::set<xtal::UnitCellCoord>>> const &
  local_orbits_sites() const {
    return m_local_orbits_sites;
  }

 private:
  /// Orbits to calculate
  std::set<int> m_orbits_to_calculate;

  /// \brief Combine orbits?
  ///
  /// If true, calculate the number of each component for
  /// the sites in the union of the orbits in `_orbits_to_calculate`. If
  /// false, calculate the number of each component for the sites in each
  /// orbit in
  /// `_orbits_to_calculate` for each orbit. If true, the resulting value
  /// will be a matrix with a single column. If false, the value will be a
  /// matrix with a column for each orbit.
  bool m_combine_orbits;

  /// Prim neighbor list
  std::shared_ptr<clexulator::PrimNeighborList> m_prim_nlist;

  /// Supercell neighbor list
  std::shared_ptr<clexulator::SuperNeighborList> m_supercell_nlist;

  /// Configuration to use
  clexulator::ConfigDoFValues const *m_dof_values;

  /// Converter from occupation index to component index, by sublattice
  /// component_index = converter[sublattice_index][occ_index]
  std::vector<std::vector<Index>> m_occ_index_to_component_index_converter;

  /// Store {neighbor_index, sublattice_index} for each site in each orbit:
  /// m_local_orbits_neighbor_indices[equivalent_index][orbit_index] ->
  ///     std::set<std::pair<int, int>>
  /// If m_combine_orbits is true, all orbits are combined into `orbit_index` 0
  std::vector<std::vector<std::set<std::pair<int, int>>>>
      m_local_orbits_neighbor_indices;

  /// Store the unique set of sites in each orbit:
  /// m_local_orbits_neighbor_indices[equivalent_index][orbit_index] ->
  ///     std::set<xta::UnitCellCoord>
  /// If m_combine_orbits is true, all orbits are combined into `orbit_index` 0
  std::vector<std::vector<std::set<xtal::UnitCellCoord>>> m_local_orbits_sites;

  /// Holds calculated composition as number of each component by orbit:
  ///
  /// If m_combine_orbits is false:
  ///     n = m_num_each_component_by_orbit(
  ///             component_index,
  ///             orbits_to_calculate_index)
  ///
  /// If m_combine_orbits is true:
  ///     n = m_num_each_component_by_orbit(component_index, 0)
  Eigen::MatrixXi m_num_each_component_by_orbit;

  /// Indices of m_local_orbits_neighbor_indices / m_local_orbits_sites to use
  /// whether combine_orbits is true or false
  std::set<int> m_unified_orbits_to_calculate;
};

}  // namespace clexmonte
}  // namespace CASM

#endif