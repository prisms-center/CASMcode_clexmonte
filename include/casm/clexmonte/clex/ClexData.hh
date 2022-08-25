#ifndef CASM_clexmonte_clex_ClexData
#define CASM_clexmonte_clex_ClexData

#include "casm/clexmonte/definitions.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/SparseCoefficients.hh"

namespace CASM {
namespace clexmonte {

/// \brief Holds data necessary for cluster expansion calculations
struct ClexData {
  ClexData(
      std::shared_ptr<clexulator::PrimNeighborList> const &_prim_neighbor_list,
      std::shared_ptr<clexulator::Clexulator> const &_clexulator,
      clexulator::SparseCoefficients const &_eci)
      : prim_neighbor_list(_prim_neighbor_list),
        clexulator(_clexulator),
        eci(_eci) {}

  /// \brief Get the neighborhood of sites where a change in DoF value requires
  ///     an update of the correlation values associated with the origin unit
  ///     cell
  std::set<xtal::UnitCellCoord> required_update_neighborhood() const {
    unsigned int const *begin = eci.index.data();
    unsigned int const *end = begin + eci.index.size();
    return clexulator->site_neighborhood(begin, end);
  }

  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;
  std::shared_ptr<clexulator::Clexulator> clexulator;
  clexulator::SparseCoefficients eci;
};

/// \brief Holds data necessary for local cluster expansion calculations
struct LocalClexData {
  LocalClexData(
      std::shared_ptr<clexulator::PrimNeighborList> const &_prim_neighbor_list,
      std::shared_ptr<std::vector<clexulator::Clexulator>> const &_clexulator,
      clexulator::SparseCoefficients const &_eci)
      : prim_neighbor_list(_prim_neighbor_list),
        clexulator(_clexulator),
        eci(_eci) {}

  /// \brief Get the neighborhood of sites where a change in DoF value requires
  ///     an update of the correlation values associated with the origin unit
  ///     cell
  std::set<xtal::UnitCellCoord> required_update_neighborhood(
      Index equivalent_index) const {
    unsigned int const *begin = eci.index.data();
    unsigned int const *end = begin + eci.index.size();
    return clexulator->at(equivalent_index).site_neighborhood(begin, end);
  }

  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;
  std::shared_ptr<std::vector<clexulator::Clexulator>> clexulator;
  clexulator::SparseCoefficients eci;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
