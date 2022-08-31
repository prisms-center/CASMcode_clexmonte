#ifndef CASM_clexmonte_events_event_data
#define CASM_clexmonte_events_event_data

#include <string>
#include <vector>

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"
#include "casm/monte/events/OccEvent.hh"

namespace CASM {
namespace clexmonte {

/// \brief Data calculated for a single event in a single state
struct EventState {
  bool is_allowed;      ///< Is allowed given current configuration
  bool is_normal;       ///< Is "normal" (dEa > 0.0) && (dEa > dEf)
  double dE_final;      ///< Final state energy, relative to initial state
  double Ekra;          ///< KRA energy
  double dE_activated;  ///< Activation energy, relative to initial state
  double freq;          ///< Attempt frequency
  double rate;          ///< Occurance rate
};

/// \brief Data particular to a single translationally distinct event
struct EventData {
  /// \brief Unit cell linear index
  Index unitcell_index;

  /// \brief Used to apply event and track occupants in monte::OccLocation
  monte::OccEvent event;
};

/// \brief Data common to all translationally equivalent events
struct PrimEventData {
  /// \brief Event type name
  std::string event_type_name;

  /// \brief Equivalent event index
  Index equivalent_index;

  /// \brief Is forward trajectory (else reverse)
  bool is_forward;

  /// \brief Linear index for this prim event
  Index prim_event_index;

  /// \brief Defines event
  occ_events::OccEvent event;

  /// \brief Event sites, relative to origin unit cell
  std::vector<xtal::UnitCellCoord> sites;

  /// \brief Initial site occupation
  std::vector<int> occ_init;

  /// \brief Final site occupation
  std::vector<int> occ_final;
};

/// \brief Data structure specifying information about the impact of a possible
///     event. Used to construct impact tables.
struct EventImpactInfo {
  /// \brief Sites whose DoF are modified when this event occurs
  std::vector<xtal::UnitCellCoord> phenomenal_sites;

  /// \brief The set of sites for which a change in DoF results in a
  ///     change in the event propensity
  std::set<xtal::UnitCellCoord> required_update_neighborhood;
};

/// \brief Identifies an event via translation from the `origin` unit cell
///
/// This data structure is used to build a "relative impact table" listing which
/// events are impacted by the occurance of each possible event in the origin
/// unit cell.
struct RelativeEventID {
  /// \brief Index specifying a possible event in the `origin` unit cell
  Index prim_event_index;

  /// \brief Translation of the event from the origin unit cell
  xtal::UnitCell translation;
};

bool operator<(RelativeEventID const &lhs, RelativeEventID const &rhs);

/// \brief Identifies an event via linear unit cell index in some supercell
///
/// Thie unitcell index and prim event index can be used to lookup the correct
/// local clexulator and neighbor list information for evaluating local
/// correlations and updating global correlations.
struct EventID {
  /// \brief Index specifying a possible event in the `origin` unit cell
  Index prim_event_index;

  /// \brief Linear unit cell index into a supercell, as determined by
  /// xtal::UnitCellIndexConverter
  Index unitcell_index;
};

bool operator<(EventID const &lhs, EventID const &rhs);

// -- Inline definitions --

inline bool operator<(RelativeEventID const &lhs, RelativeEventID const &rhs) {
  if (lhs.translation < rhs.translation) {
    return true;
  }
  if (rhs.translation < lhs.translation) {
    return false;
  }
  return lhs.prim_event_index < rhs.prim_event_index;
}

inline bool operator<(EventID const &lhs, EventID const &rhs) {
  if (lhs.unitcell_index < rhs.unitcell_index) {
    return true;
  }
  if (lhs.unitcell_index > rhs.unitcell_index) {
    return false;
  }
  return lhs.prim_event_index < rhs.prim_event_index;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
