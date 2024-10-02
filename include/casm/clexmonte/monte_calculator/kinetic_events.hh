#ifndef CASM_clexmonte_monte_calculator_kinetic_events
#define CASM_clexmonte_monte_calculator_kinetic_events

#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

class LocalOrbitCompositionCalculator {
 public:
  LocalOrbitCompositionCalculator(std::shared_ptr<system_type> _system,
                                  std::string _event_type_name,
                                  std::set<int> _orbits_to_calculate);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state);

  /// \brief Calculate the composition by orbit around an event
  Eigen::MatrixXi const &calculate_num_each_component(
      Eigen::VectorXi const &occupation, Index unitcell_index,
      Index equivalent_index);

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// Event type name
  std::string m_event_type_name;

  /// Orbits to calculate
  std::set<int> m_orbits_to_calculate;

  /// State to use
  state_type const *m_state;

  /// Converter from occupation index to component index, by sublattice
  /// component_index = converter[sublattice_index][occ_index]
  std::vector<std::vector<Index>> m_occ_index_to_component_index_converter;

  /// Holds calculated composition as number of each component by orbit:
  ///
  ///     n = m_num_each_component_by_orbit(
  ///             component_index,
  ///             orbits_to_calculate_index)
  ///
  Eigen::MatrixXi m_num_each_component_by_orbit;

  /// Supercell neighbor list
  std::shared_ptr<clexulator::SuperNeighborList> m_supercell_nlist;

  /// Store {neighbor_index, sublattice_index} for each site in each orbit:
  /// m_local_orbits_neighbor_indices[equivalent_index][orbit_index] ->
  ///     std::set<std::pair<int, int>>
  std::vector<std::vector<std::set<std::pair<int, int>>>>
      m_local_orbits_neighbor_indices;
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
  /// \brief Constructor
  EventStateCalculator(std::shared_ptr<system_type> _system,
                       std::string _event_type_name);

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state);

  /// \brief Pointer to current state
  state_type const *state() const;

  /// \brief Current state's reciprocal temperature
  double beta() const;

  /// \brief Calculate the state of an event
  void calculate_event_state(EventState &state, EventData const &event_data,
                             PrimEventData const &prim_event_data) const;

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients() const {
    if (m_formation_energy_clex == nullptr) {
      throw std::runtime_error(
          "EventStateCalculator::formation_energy_coefficients: "
          "m_formation_energy_clex == nullptr");
    }
    return m_formation_energy_clex->coefficients();
  }

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
  std::shared_ptr<LocalBasisSetClusterInfo> m_cluster_info;
};

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
  Index not_normal_count;

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
class KineticEventData : public BaseMonteEventData {
 public:
  typedef std::mt19937_64 engine_type;
  typedef lotto::RejectionFreeEventSelector<EventID, CompleteEventCalculator,
                                            engine_type>
      event_selector_type;

  KineticEventData(std::shared_ptr<system_type> _system,
                   std::optional<std::vector<EventFilterGroup>> _event_filters);

  // --- BaseMonteEventData data ---

  /// The system
  // std::shared_ptr<system_type> system;
  using BaseMonteEventData::system;

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  // std::vector<clexmonte::PrimEventData> prim_event_list;
  using BaseMonteEventData::prim_event_list;

  /// Information about what sites may impact each prim event
  // std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;
  using BaseMonteEventData::prim_impact_info_list;

  // -- Data - set when `update` is called --

  /// Functions for calculating event states, one for each prim event.
  /// This is supercell-specific, even though it is one per prim event,
  /// because it depends on supercell-specific clexulators
  std::vector<EventStateCalculator> prim_event_calculators;

  /// Current supercell
  Eigen::Matrix3l transformation_matrix_to_super;

  /// Selectively allow events by unit cell
  std::vector<EventFilterGroup> event_filters;

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Calculator for KMC event selection
  std::shared_ptr<CompleteEventCalculator> event_calculator;

  /// Event selector
  std::shared_ptr<event_selector_type> event_selector;

  /// \brief Update for given state, conditions, and occupants
  void update(state_type const &state, monte::OccLocation const &occ_location,
              std::optional<std::vector<EventFilterGroup>> _event_filters,
              std::shared_ptr<engine_type> engine);

  // --- Event selection ---

  /// \brief Select an event, and optionally re-calculate event state for the
  ///     selected event
  void select_event(SelectedEvent &selected_event, bool requires_event_state);

  // --- BaseMonteEventData interface ---

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients()
      const override {
    if (prim_event_calculators.size() == 0) {
      throw std::runtime_error(
          "KineticEventData::formation_energy_coefficients: "
          "prim_event_calculators.size() == 0");
    }
    return prim_event_calculators.at(0).formation_energy_coefficients();
  }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const override {
    if (prim_event_calculators.size() == 0) {
      throw std::runtime_error(
          "KineticEventData::kra_coefficients: "
          "prim_event_calculators.size() == 0");
    }
    if (prim_event_index >= prim_event_calculators.size()) {
      throw std::runtime_error(
          "KineticEventData::kra_coefficients: "
          "prim_event_index >= prim_event_calculators.size()");
    }
    return prim_event_calculators.at(prim_event_index).freq_coefficients();
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const override {
    if (prim_event_calculators.size() == 0) {
      throw std::runtime_error(
          "KineticEventData::kra_coefficients: "
          "prim_event_calculators.size() == 0");
    }
    if (prim_event_index >= prim_event_calculators.size()) {
      throw std::runtime_error(
          "KineticEventData::kra_coefficients: "
          "prim_event_index >= prim_event_calculators.size()");
    }
    return prim_event_calculators.at(prim_event_index).kra_coefficients();
  }

  // -- Event list summary info --

  /// The size of the event list
  Index n_events() const override { return event_list.events.size(); }

  /// Return the current total event rate
  double total_rate() const override {
    if (event_selector == nullptr) {
      throw std::runtime_error(
          "KineticEventData::total_rate: Events have not been calculated");
    }
    return event_selector->total_rate();
  }

  // -- Event list iteration --

  /// Construct new internal iterator and return its index
  Index new_iterator(bool is_end) override {
    Index i = 0;
    while (m_it.find(i) != m_it.end()) {
      ++i;
    }
    if (is_end) {
      m_it.emplace(i, event_list.events.end());
    } else {
      m_it.emplace(i, event_list.events.begin());
    }
    return i;
  }

  /// Copy internal iterator and return the new iterator index
  Index copy_iterator(Index i) override {
    if (m_it.find(i) == m_it.end()) {
      throw std::runtime_error(
          "KineticEventData::copy_iterator: Iterator not found");
    }
    Index j = 0;
    while (m_it.find(j) != m_it.end()) {
      ++j;
    }
    m_it.emplace(j, m_it[i]);
    return j;
  }

  /// Erase internal iterator
  void erase_iterator(Index i) override { m_it.erase(i); }

  /// Check if two internal iterators are equal
  bool equal_iterator(Index i, Index j) override {
    auto it_i = m_it.find(i);
    auto it_j = m_it.find(j);
    if (it_i == m_it.end() || it_j == m_it.end()) {
      throw std::runtime_error(
          "KineticEventData::equal_iterator: Iterator not found");
    }
    return it_i->second == it_j->second;
  }

  /// Advance internal iterator by one event
  void advance_iterator(Index i) override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "KineticEventData::advance_iterator: Iterator not found");
    }
    if (it->second == event_list.events.end()) {
      throw std::runtime_error(
          "KineticEventData::advance_iterator: Cannot advance past end of "
          "event list");
    }
    ++it->second;
  }

  /// The event ID for the current state of the internal iterator
  EventID const &event_id(Index i) const override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "KineticEventData::event_id: Iterator not found");
    }
    return it->second->first;
  }

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the specified event. Reference is
  /// valid until the next call to this method.
  monte::OccEvent const &event_to_apply(EventID const &id) const override {
    auto it = event_list.events.find(id);
    if (it == event_list.events.end()) {
      throw std::runtime_error(
          "KineticEventData::event_to_apply: Event not found in event list");
    }
    return it->second.event;
  }

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const override {
    if (event_selector == nullptr) {
      throw std::runtime_error(
          "KineticEventData::total_rate: Events have not been calculated");
    }
    return event_selector->get_rate(id);
  }

  /// Calculate event state data
  ///
  /// Notes:
  /// - Event state is only valid for event `id` until the next call to this
  ///   method
  EventState const &event_state(EventID const &id) const override {
    auto it = event_list.events.find(id);
    if (it == event_list.events.end()) {
      throw std::runtime_error(
          "KineticEventData::event_state: Event not found in event list");
    }
    EventData const &_event_data = it->second;
    PrimEventData const &_prim_event_data =
        prim_event_list.at(id.prim_event_index);
    prim_event_calculators.at(id.prim_event_index)
        .calculate_event_state(m_event_state, _event_data, _prim_event_data);
    return m_event_state;
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &impact(EventID const &id) const override {
    auto it = event_list.impact_table.find(id);
    if (it == event_list.impact_table.end()) {
      throw std::runtime_error(
          "KineticEventData::impact: Event not found in impact table");
    }
    return it->second;
  }

 private:
  /// \brief Holds internal iterators to allow generic interface for
  /// iterating over EventID
  std::map<Index, std::map<EventID, EventData>::const_iterator> m_it;

  /// \brief Holds temporary calculated event state
  mutable EventState m_event_state;
};

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

#endif