#ifndef CASM_clexmonte_BaseMonteEventData
#define CASM_clexmonte_BaseMonteEventData

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

namespace CASM {
namespace clexmonte {

struct EventFilterGroup;
struct StateData;

class EventStateCalculator;

typedef std::function<void(std::reference_wrapper<EventState> state,
                           EventStateCalculator const &calculator)>
    CustomEventStateCalculationFunction;

/// \brief Event rate calculation for a particular KMC event
///
/// EventStateCalculator is used to separate the event calculation from the
/// event definition data in PrimEventData. All symmetrically equivalent
/// events can use the same EventStateCalculator, but a simple approach
/// is to create one for each distinct event associated with the primitive
/// cell.
class EventStateCalculator {
 public:
  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state);

  /// \brief Set custom event state calculation function
  void set_custom_event_state_calculation(
      CustomEventStateCalculationFunction f);

  /// \brief Clear custom event state calculation function
  void clear_custom_event_state_calculation();

  /// \brief Pointer to current state
  state_type const *state() const { return m_state; }

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, Index unitcell_index,
                             std::vector<Index> const &linear_site_index,
                             PrimEventData const &prim_event_data) const;

  /// \brief Return the event type name
  std::string const &event_type_name() const { return m_event_type_name; }

  /// \brief Current state's temperature
  double temperature() const { return *m_temperature; }

  /// \brief Current state's reciprocal temperature
  double beta() const { return 1.0 / (CASM::KB * *this->m_temperature); }

  /// \brief Get the unitcell index for the event currently being calculated
  /// (valid for custom event state calculations only)
  Index curr_unitcell_index() const { return m_unitcell_index; }

  /// \brief Get the linear site indices for the event currently being
  /// calculated (valid for custom event state calculations only)
  std::vector<Index> const &curr_linear_site_index() const {
    return *m_linear_site_index;
  }

  /// \brief Get the prim event data for the event currently being calculated
  /// (valid for custom event state calculations only)
  PrimEventData const &curr_prim_event_data() const {
    return *m_prim_event_data;
  }

  /// Set default event state (short cut for custom event state calculations)
  void set_default_event_state(EventState &state) const {
    _default_event_state_calculation(state, this->curr_unitcell_index(),
                                     this->curr_linear_site_index(),
                                     this->curr_prim_event_data());
  }

  /// Get the formation energy cluster expansion
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex() const {
    return m_formation_energy_clex;
  }

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients() const {
    if (m_formation_energy_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::formation_energy_coefficients: "
          "m_formation_energy_clex == nullptr");
    }
    return m_formation_energy_clex->coefficients();
  }

  /// Get the event multi-local cluster expansion
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> event_clex() const {
    return m_event_clex;
  }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the KRA value
  Index kra_index() const { return m_kra_index; }

  /// The index of the event multi-local cluster expansion output that
  /// corresponds to the attempt frequency value
  Index freq_index() const { return m_freq_index; }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::freq_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_freq_index];
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients() const {
    if (m_event_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::kra_coefficients: m_event_clex == nullptr");
    }
    return m_event_clex->coefficients()[m_kra_index];
  }

 private:
  /// \brief Calculate the state of an event
  void _default_event_state_calculation(
      EventState &state, Index unitcell_index,
      std::vector<Index> const &linear_site_index,
      PrimEventData const &prim_event_data) const;

  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// State to use
  state_type const *m_state;

  /// Current state's temperature
  double const *m_temperature;

  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
  mutable Eigen::VectorXd m_event_values;
  Index m_kra_index;
  Index m_freq_index;

  /// If true, use custom event state calculation function
  bool m_custom_event_state_calculation;

  /// Custom event state calculation function
  CustomEventStateCalculationFunction m_custom_event_state_calculation_f;

  /// Current event's unitcell index
  mutable Index m_unitcell_index;

  /// Current event's linear site indices
  mutable std::vector<Index> const *m_linear_site_index;

  /// Current event's prim event data
  mutable PrimEventData const *m_prim_event_data;
};

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

  /// Custom event state calculation functions
  std::map<std::string, CustomEventStateCalculationFunction>
      custom_event_state_calculation_f;

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

  /// \brief Set a custom event state calculation function
  void set_custom_event_state_calculation(
      std::string const &event_type_name,
      CustomEventStateCalculationFunction f) {
    custom_event_state_calculation_f[event_type_name] = f;
  }

  /// \brief Erase a custom event state calculation function
  void erase_custom_event_state_calculation(
      std::string const &event_type_name,
      CustomEventStateCalculationFunction f) {
    custom_event_state_calculation_f.erase(event_type_name);
  }

  virtual void update(
      std::shared_ptr<StateData> _state_data,
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
