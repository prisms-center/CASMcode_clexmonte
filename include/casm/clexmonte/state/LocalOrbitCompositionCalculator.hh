#ifndef CASM_clexmonte_LocalOrbitCompositionCalculator
#define CASM_clexmonte_LocalOrbitCompositionCalculator

#include "casm/clexulator/NeighborList.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/LinearIndexConverter.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
}

namespace clexmonte {
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
      std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
      xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
      composition::CompositionCalculator const &_composition_calculator,
      clexulator::ConfigDoFValues const *_dof_values = nullptr);

  /// \brief Reset pointer to configuration currently being calculated
  void set(clexulator::ConfigDoFValues const *dof_values);

  // {prim_event_index} -> \vec{n}
  // {equivalent_index} -> \vec{n}
  // {event_type} -> \vec{n}

  /// \brief Calculate the composition by orbit around an event
  Eigen::MatrixXi const &value(Index unitcell_index, Index equivalent_index);

 private:
  /// Orbits to calculate
  std::set<int> m_orbits_to_calculate;

  /// \brief Combine orbits?
  ///
  /// If true, calculate the number of each component for
  /// the sites in the union of the orbits in `_orbits_to_calculate`. If false,
  /// calculate the number of each component for the sites in each orbit in
  /// `_orbits_to_calculate` for each orbit. If true, the resulting value will
  /// be a matrix with a single column. If false, the value will be a matrix
  /// with a column for each orbit.
  bool m_combine_orbits;

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
  std::vector<std::vector<std::set<std::pair<int, int>>>>
      m_local_orbits_neighbor_indices;

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
};

}  // namespace clexmonte
}  // namespace CASM

#endif