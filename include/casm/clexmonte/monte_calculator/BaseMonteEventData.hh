#ifndef CASM_clexmonte_BaseMonteEventData
#define CASM_clexmonte_BaseMonteEventData

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/sampling/SelectedEventData.hh"

namespace CASM {
namespace clexmonte {

struct EventFilterGroup;

/// \brief Base class to provide access to event data for a Monte Carlo
/// simulation
class BaseMonteEventData {
 public:
  typedef std::mt19937_64 engine_type;
  typedef monte::KMCData<config_type, statistics_type, engine_type>
      kmc_data_type;
  typedef clexmonte::run_manager_type<engine_type> run_manager_type;

  BaseMonteEventData() = default;
  virtual ~BaseMonteEventData() = default;

  /// The system
  std::shared_ptr<system_type> system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  // -- System data --

  /// Get the formation energy coefficients
  virtual clexulator::SparseCoefficients const &formation_energy_coefficients()
      const = 0;

  /// Get the attempt frequency coefficients for a specific event
  virtual clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const = 0;

  /// Get the KRA coefficients for a specific event
  virtual clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const = 0;

  // -- Update and run --

  virtual void update(
      state_type const &state, monte::OccLocation const &occ_location,
      std::optional<std::vector<EventFilterGroup>> _event_filters,
      std::shared_ptr<engine_type> engine) = 0;

  virtual void run(state_type &state, monte::OccLocation &occ_location,
                   kmc_data_type &kmc_data, SelectedEvent &selected_event,
                   std::optional<monte::SelectedEventDataCollector> &collector,
                   run_manager_type &run_manager,
                   std::shared_ptr<occ_events::OccSystem> event_system) = 0;

  // -- Select Event --

  /// Select an event to apply
  virtual void select_event(SelectedEvent &selected_event,
                            bool requires_event_state) = 0;

  /// Return number of events calculated with no barrier, by type
  virtual std::map<std::string, Index> const &n_not_normal() const = 0;

  // -- Event list summary info --

  /// The size of the event list
  virtual Index n_events() const = 0;

  /// Return the current total event rate
  virtual double total_rate() const = 0;

  // -- Event list iteration --

  /// Construct new internal iterator and return its index
  virtual Index new_iterator(bool is_end) = 0;

  /// Copy internal iterator and return the new iterator index
  virtual Index copy_iterator(Index i) = 0;

  /// Erase internal iterator
  virtual void erase_iterator(Index i) = 0;

  /// Check if two internal iterators are equal
  virtual bool equal_iterator(Index i, Index j) = 0;

  /// Advance internal iterator by one event
  virtual void advance_iterator(Index i) = 0;

  /// The event ID for the current state of the internal iterator
  virtual EventID const &event_id(Index i) const = 0;

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the specified event. Reference is
  /// valid until the next call to this method.
  virtual monte::OccEvent const &event_to_apply(EventID const &id) const = 0;

  /// Return the current rate for a specific event
  virtual double event_rate(EventID const &id) const = 0;

  /// Calculate event state data. Reference is valid until the next call to this
  /// method.
  virtual EventState const &event_state(EventID const &id) const = 0;

  /// The events that must be updated if the specified event occurs. Reference
  /// is valid until the next call to this method.
  virtual std::vector<EventID> const &impact(EventID const &id) const = 0;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
