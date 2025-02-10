#ifndef CASM_clexmonte_kinetic_events
#define CASM_clexmonte_kinetic_events

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/LocalClusterExpansion.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

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

/// \brief Event rate calculation for a particular KMC event
///
/// EventStateCalculator is used to separate the event calculation from the
/// event definition data in PrimEventData. All symmetrically equivalent
/// events can use the same EventStateCalculator, but a simple approach
/// is to create one for each distinct event associated with the primitive
/// cell.
class EventStateCalculator {
 public:
  // /// \brief Constructor
  // EventStateCalculator(
  //     std::shared_ptr<Conditions> _conditions,
  //     std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex,
  //     std::shared_ptr<clexulator::MultiLocalClusterExpansion> _event_clex,
  //     std::string _name, std::map<std::string, Index> _glossary);

  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state, std::shared_ptr<Conditions> conditions);

  /// \brief Pointer to current state
  state_type const *state() const;

  /// \brief Pointer to current conditions
  std::shared_ptr<Conditions> const &conditions() const;

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, EventData const &event_data,
                             PrimEventData const &prim_event_data) const;

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// State to use
  state_type const *m_state;

  /// Conditions
  std::shared_ptr<Conditions> m_conditions;

  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
  Index m_kra_index;
  Index m_freq_index;
};

/// \brief Construct a vector EventStateCalculator, one per event in a
///     vector of PrimEventData
std::vector<EventStateCalculator> make_prim_event_calculators(
    std::shared_ptr<system_type> system, state_type const &state,
    std::vector<PrimEventData> const &prim_event_list,
    std::shared_ptr<Conditions> conditions);

/// \brief CompleteEventCalculator is an event calculator with the required
/// interface for the
///     classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct CompleteEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Prim event calculators - order must match prim_event_list
  std::vector<EventStateCalculator> const &prim_event_calculators;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Write to warn about non-normal events
  Log &event_log;

  // Note: to keep all event state calculations, comment out this:
  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Count not-normal events
  Index abnormal_count;

  CompleteEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<EventStateCalculator> const &_prim_event_calculators,
      std::map<EventID, EventData> const &_event_list,
      Log &_event_log = CASM::err_log());

  /// \brief Get CASM::monte::OccEvent corresponding to given event ID
  double calculate_rate(EventID const &id);
};

/// \brief Data for kinetic Monte Carlo events
///
/// Includes:
/// - prim event list
/// - prim impact info
/// - event state calculators: one per prim event, given a state pointer and
///   can then calculate event energies, attempt frequency, and rate for the
///   for the current state on request
/// - complete event list
/// - CompleteEventCalculator: uses event state calculator and complete event
///   list to calculate a rate given an event ID
struct KineticEventData {
  KineticEventData(std::shared_ptr<system_type> _system);

  /// \brief Update for given state, conditions, and occupants
  void update(state_type const &state, std::shared_ptr<Conditions> conditions,
              monte::OccLocation const &occ_location,
              std::vector<EventFilterGroup> const &event_filters);

  /// The system
  std::shared_ptr<system_type> system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Functions for calculating event states, one for each prim event.
  /// This is supercell-specific, even though it is one per prim event,
  /// because it depends on supercell-specific clexulators
  std::vector<EventStateCalculator> prim_event_calculators;

  /// Calculator for KMC event selection
  std::shared_ptr<CompleteEventCalculator> event_calculator;
};

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
