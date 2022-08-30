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
  bool is_allowed;  ///< Is allowed given current configuration
  bool is_normal;   ///< Is "normal" (dEa > 0.0) && (dEa > dEf)
  double dEf;       ///< Final state energy, relative to initial state
  double Ekra;      ///< KRA energy
  double dEa;       ///< Activation energy, relative to initial state
  double freq;      ///< Attempt frequency
  double rate;      ///< Occurance rate
};

/// \brief Data particular to a single translationally distinct event
struct EventData {
  /// \brief Unit cell linear index
  Index unitcell_index;

  /// \brief Used to apply event and track occupants in monte::OccLocation
  monte::OccEvent event;

  /// \brief Calculated event properties in current state
  EventState event_state;
};

/// \brief Data common to all translationally equivalent events
struct PrimEventData {
  /// \brief Event type name
  std::string event_type_name;

  /// \brief Equivalent event index
  Index equivalent_index;

  /// \brief Is forward trajectory (else reverse)
  bool is_forward;

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

}  // namespace clexmonte
}  // namespace CASM

#endif
