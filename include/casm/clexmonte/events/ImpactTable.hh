#ifndef CASM_clexmonte_events_ImpactTable
#define CASM_clexmonte_events_ImpactTable

#include <set>
#include <vector>

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace clexmonte {

/// \brief Implements an event impact table, storing only relative interations
/// explicitly
///
///
/// Specifies which events are impacted (and therefore must have their
/// propensities updated) by the occurance of another event.
///
/// RelativeEventImpactTable generates impact vectors for a particular event
/// occurance upon request by applying translations and within
/// supercell operations. Therefore, compared to SupercellEventImpactTable it
/// has lower memory requirements but will be slower.
struct RelativeEventImpactTable {
  RelativeEventImpactTable(
      std::vector<EventImpactInfo> const &prim_event_list,
      xtal::UnitCellIndexConverter const &unitcell_converter);

  std::vector<EventID> const &operator()(EventID const &event_id) const;

 private:
  std::vector<std::vector<RelativeEventID>> m_impact_table;
  xtal::UnitCellIndexConverter m_unitcell_converter;
  mutable std::vector<EventID> m_result;
};

/// \brief Implements an event impact table, storing only relative interations
/// explicitly, and using the neighbor list to determine unitcell indices
struct NeighborlistEventImpactTable {
  NeighborlistEventImpactTable(
      std::vector<EventImpactInfo> const &prim_event_list,
      xtal::UnitCellIndexConverter const &unitcell_converter,
      std::shared_ptr<clexulator::PrimNeighborList> prim_nlist,
      Eigen::Matrix3l const &transformation_matrix_to_super,
      std::shared_ptr<clexulator::SuperNeighborList> supercell_nlist);

  std::vector<EventID> const &operator()(EventID const &event_id) const;

 private:
  /// \brief Neighbor list indices used to set the unitcell_index of impacted
  /// events in `results`
  std::vector<std::vector<Index>> m_unitcell_nlist_index;
  std::shared_ptr<clexulator::SuperNeighborList> m_supercell_nlist;
  mutable std::vector<std::vector<EventID>> m_result;
};

/// \brief Implement an event impact table, storing all interations in a
/// supercell explicitly
///
/// Specifies which events are impacted (and therefore must have their
/// propensities updated) by the occurance of another event.
///
/// SupercellEventImpactTable generates and stores all impact vectors at
/// construction. Therefore, compared to RelativeEventImpactTable it
/// has higher memory requirements but will be faster.
struct SupercellEventImpactTable {
  SupercellEventImpactTable(
      std::vector<EventImpactInfo> const &prim_event_list,
      xtal::UnitCellIndexConverter const &unitcell_converter);

  std::vector<EventID> const &operator()(EventID const &event_id) const;

 private:
  Index m_n_prim_events;
  std::vector<std::vector<EventID>> m_impact_table;
};

/// \brief Return an impact table for events in the origin unit cell
std::vector<std::vector<RelativeEventID>> make_relative_impact_table(
    std::vector<EventImpactInfo> const &prim_event_list);

// -- Inline definitions --

inline std::vector<EventID> const &RelativeEventImpactTable::operator()(
    EventID const &event_id) const {
  std::vector<RelativeEventID> const &relative_impact =
      m_impact_table[event_id.prim_event_index];

  m_result.resize(relative_impact.size());
  for (Index i = 0; i < relative_impact.size(); ++i) {
    RelativeEventID const &relative_event_id = relative_impact[i];
    m_result[i].prim_event_index = relative_event_id.prim_event_index;
    m_result[i].unitcell_index =
        m_unitcell_converter(m_unitcell_converter(event_id.unitcell_index) +
                             relative_event_id.translation);
  }
  return m_result;
}

inline std::vector<EventID> const &NeighborlistEventImpactTable::operator()(
    EventID const &event_id) const {
  std::vector<Index> const &nlist_index =
      m_unitcell_nlist_index[event_id.prim_event_index];
  std::vector<Index> const &unitcells =
      m_supercell_nlist->unitcells(event_id.unitcell_index);
  std::vector<EventID> &result = m_result[event_id.prim_event_index];
  for (Index i = 0; i < nlist_index.size(); ++i) {
    result[i].unitcell_index = unitcells[nlist_index[i]];
  }
  return result;
}

inline std::vector<EventID> const &SupercellEventImpactTable::operator()(
    EventID const &event_id) const {
  Index linear_index =
      event_id.unitcell_index * m_n_prim_events + event_id.prim_event_index;
  return m_impact_table[linear_index];
}

}  // namespace clexmonte
}  // namespace CASM

#endif
