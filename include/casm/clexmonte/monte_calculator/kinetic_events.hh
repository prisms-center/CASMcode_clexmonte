#ifndef CASM_clexmonte_monte_calculator_kinetic_events
#define CASM_clexmonte_monte_calculator_kinetic_events

#include "casm/clexmonte/events/AllowedEventList.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/lotto.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

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
  void calculate_event_state(EventState &state, Index unitcell_index,
                             std::vector<Index> const &linear_site_index,
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

  /// \brief Count not-normal events (key == event_type_name; value == count)
  std::map<std::string, Index> n_not_normal;

  CompleteEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<EventStateCalculator> const &_prim_event_calculators,
      std::map<EventID, EventData> const &_event_list,
      Log &_event_log = CASM::err_log());

  /// \brief Update `event_state` for event `id` in the current state and
  /// return the event rate
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
template <bool DebugMode>
class CompleteKineticEventData : public BaseMonteEventData {
 public:
  typedef std::mt19937_64 engine_type;
  typedef lotto::RejectionFreeEventSelector<EventID, CompleteEventCalculator,
                                            engine_type>
      event_selector_type;

  CompleteKineticEventData(
      std::shared_ptr<system_type> _system,
      std::optional<std::vector<EventFilterGroup>> _event_filters,
      bool _allow_events_with_no_barrier);

  // -- Options --

  /// \brief If false (default), events with no barrier are not allowed and
  /// an exception is thrown by `select_event` if one was encountered when
  /// calculating rates in the preceding step. If true, events with no barrier
  /// are allowed with rate 1.0 and warning messages are written.
  bool allow_events_with_no_barrier = false;

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
              std::shared_ptr<engine_type> engine) override;

  /// \brief Run the KMC simulation
  void run(state_type &state, monte::OccLocation &occ_location,
           kmc_data_type &kmc_data, SelectedEvent &selected_event,
           std::optional<monte::SelectedEventDataCollector> &collector,
           run_manager_type &run_manager,
           std::shared_ptr<occ_events::OccSystem> event_system) override;

  // -- Validators --

  EventStateCalculator const &_prim_event_calculator(
      Index prim_event_index) const {
    if (prim_event_calculators.size() == 0) {
      throw std::runtime_error(
          "Error in CompleteKineticEventData: "
          "prim_event_calculators.size() == 0");
    }
    if (prim_event_index >= prim_event_calculators.size()) {
      throw std::runtime_error(
          "CompleteKineticEventData::kra_coefficients: "
          "prim_event_index (=" +
          std::to_string(prim_event_index) +
          ") >= prim_event_calculators.size()");
    }
    return prim_event_calculators.at(prim_event_index);
  }

  std::pair<const EventID, EventData> const &_event(EventID const &id) const {
    auto it = event_list.events.find(id);
    if (it == event_list.events.end()) {
      std::stringstream ss;
      ss << "Error in CompleteKineticEventData: Event (prim_event_index="
         << id.prim_event_index << ", unitcell_index=" << id.unitcell_index
         << ") not found in event list";
      throw std::runtime_error(ss.str());
    }
    return *it;
  }

  CompleteEventCalculator &_event_calculator() const {
    if (!event_calculator) {
      throw std::runtime_error(
          "Error in CompleteKineticEventData: Event calculator not set");
    }
    return *event_calculator;
  }

  event_selector_type &_event_selector() const {
    if (!event_selector) {
      throw std::runtime_error(
          "Error in CompleteKineticEventData: Event selector not set");
    }
    return *event_selector;
  }

  // --- BaseMonteEventData interface ---

  // -- System data --

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients()
      const override {
    return _prim_event_calculator(0).formation_energy_coefficients();
  }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const override {
    return _prim_event_calculator(prim_event_index).freq_coefficients();
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const override {
    return _prim_event_calculator(prim_event_index).kra_coefficients();
  }

  // --- Event selection ---

  /// \brief Select an event, and optionally re-calculate event state for the
  ///     selected event
  void select_event(SelectedEvent &selected_event,
                    bool requires_event_state) override;

  /// Return number of events calculated with no barrier, by type
  std::map<std::string, Index> const &n_not_normal() const override {
    return _event_calculator().n_not_normal;
  }

  // -- Event list summary info --

  /// The size of the event list
  Index n_events() const override { return event_list.events.size(); }

  /// Return the current total event rate
  double total_rate() const override { return _event_selector().total_rate(); }

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
          "CompleteKineticEventData::copy_iterator: Iterator not found");
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
          "CompleteKineticEventData::equal_iterator: Iterator not found");
    }
    return it_i->second == it_j->second;
  }

  /// Advance internal iterator by one event
  void advance_iterator(Index i) override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "CompleteKineticEventData::advance_iterator: Iterator not found");
    }
    if (it->second == event_list.events.end()) {
      throw std::runtime_error(
          "CompleteKineticEventData::advance_iterator: "
          "Cannot advance past end of event list");
    }
    ++it->second;
  }

  /// The event ID for the current state of the internal iterator
  EventID const &event_id(Index i) const override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "CompleteKineticEventData::event_id: Iterator not found");
    }
    return it->second->first;
  }

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the specified event. Reference is
  /// valid until the next call to this method.
  monte::OccEvent const &event_to_apply(EventID const &id) const override {
    return _event(id).second.event;
  }

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const override {
    return _event_selector().get_rate(id);
  }

  /// Calculate event state data
  ///
  /// Notes:
  /// - Event state is only valid for event `id` until the next call to this
  ///   method
  EventState const &event_state(EventID const &id) const override {
    auto const &event = _event(id);
    EventData const &_event_data = event.second;
    PrimEventData const &_prim_event_data =
        prim_event_list.at(id.prim_event_index);
    prim_event_calculators.at(id.prim_event_index)
        .calculate_event_state(m_event_state, _event_data.unitcell_index,
                               _event_data.event.linear_site_index,
                               _prim_event_data);
    return m_event_state;
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &impact(EventID const &id) const override {
    auto it = event_list.impact_table.find(id);
    if (it == event_list.impact_table.end()) {
      throw std::runtime_error(
          "CompleteKineticEventData::impact: Event not found in impact table");
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

/// \brief AllowedEventCalculator is an event calculator with the required
/// interface for the
///     classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct AllowedEventCalculator {
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Prim event calculators - order must match prim_event_list
  std::vector<EventStateCalculator> const &prim_event_calculators;

  /// \brief Allowed event list
  AllowedEventList &event_list;

  /// \brief Write to warn about non-normal events
  Log &event_log;

  // Note: to keep all event state calculations, comment out this:
  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Count not-normal events (key == event_type_name; value == count)
  std::map<std::string, Index> n_not_normal;

  /// \brief Event site linear indices
  ///
  /// Used temporarily to calculate event state
  std::vector<Index> linear_site_index;

  /// \brief Event data
  ///
  /// Used temporarily to set the monte::OccEvent used to apply a selected
  /// event
  EventData event_data;

  AllowedEventCalculator(
      std::vector<PrimEventData> const &_prim_event_list,
      std::vector<EventStateCalculator> const &_prim_event_calculators,
      AllowedEventList &_event_list, Log &_event_log = CASM::err_log());

  /// \brief Update `event_state` for event `event_index` in the current state
  /// and return the event rate; if the event is no longer allowed, free the
  /// event.
  double calculate_rate(Index event_index);

  /// \brief Update `event_state` for any event `event_id` in the current state
  /// and return the event rate
  double calculate_rate(EventID const &event_id);

  /// \brief Set `event_data` for event `event_index`, returning a reference
  /// which is valid until the next call to this method
  EventData const &set_event_data(Index event_index);

  /// \brief Set `event_data` for any event `event_id`, returning a reference
  /// which is valid until the next call to this method
  EventData const &set_event_data(EventID const &event_id);
};

enum class kinetic_event_selector_type {
  vector_sum_tree,
  sum_tree,
  direct_sum,
};

typedef lotto::VectorRejectionFreeEventSelector<Index, AllowedEventCalculator,
                                                std::mt19937_64,
                                                GetImpactFromAllowedEventList>
    vector_sum_tree_event_selector_type;

typedef lotto::RejectionFreeEventSelector<Index, AllowedEventCalculator,
                                          std::mt19937_64,
                                          GetImpactFromAllowedEventList>
    sum_tree_event_selector_type;

typedef lotto::DirectSumRejectionFreeEventSelector<
    Index, AllowedEventCalculator, std::mt19937_64,
    GetImpactFromAllowedEventList>
    direct_sum_event_selector_type;

/// \brief Data for kinetic Monte Carlo events
///
/// Includes:
/// - prim event list
/// - prim impact info
/// - event state calculators: one per prim event, given a state pointer and
///   can then calculate event energies, attempt frequency, and rate for the
///   for the current state on request
/// - allowed event list
/// - AllowedEventCalculator: uses event state calculator and allowed event
///   list to calculate a rate given an event index
template <typename EventSelectorType, bool DebugMode>
class AllowedKineticEventData : public BaseMonteEventData {
 public:
  typedef std::mt19937_64 engine_type;
  typedef EventSelectorType event_selector_type;

  AllowedKineticEventData(std::shared_ptr<system_type> _system,
                          bool _allow_events_with_no_barrier = false,
                          bool _use_map_index = true,
                          bool _use_neighborlist_impact_table = true,
                          bool _assign_allowed_events_only = true);

  // -- Options --

  /// \brief If false (default), events with no barrier are not allowed and
  /// an exception is thrown by `select_event` if one was encountered when
  /// calculating rates in the preceding step. If true, events with no barrier
  /// are allowed with rate 1.0 and warning messages are written.
  const bool allow_events_with_no_barrier;

  /// \brief If true (default), use the map index for the AllowedEventMap; If
  /// false, use the vector index
  ///
  /// The map index is lower memory, but may be slower; The vector index is
  /// higher memory, but may be faster.
  const bool use_map_index;

  /// \brief If true (default), use the neighborlist impact table; else use the
  /// relative impact table
  const bool use_neighborlist_impact_table;

  /// \brief If true (default) check if potentially impacted events are allowed
  ///     and only assign them to the event list if they are. Otherwise,
  ///     assign all potentially impacted events to the event list (whether they
  ///     are allowed will still be checked during the rate calculation).
  ///
  const bool assign_allowed_events_only;

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

  /// \brief Random number generator
  ///
  /// This is constructed at `update` and stored to allow re-building the
  /// the event selector with the same random number generator
  std::shared_ptr<lotto::RandomGenerator> random_generator;

  /// Functions for calculating event states, one for each prim event.
  /// This is supercell-specific, even though it is one per prim event,
  /// because it depends on supercell-specific clexulators
  std::vector<EventStateCalculator> prim_event_calculators;

  /// All supercell events, and which events must be updated
  /// when one occurs
  std::shared_ptr<clexmonte::AllowedEventList> event_list;

  /// Calculator for KMC event selection
  std::shared_ptr<AllowedEventCalculator> event_calculator;

  // -- Event selector options --

  /// \brief Event selector
  std::shared_ptr<event_selector_type> event_selector;

  /// \brief Update for given state, conditions, and occupants
  void update(state_type const &state, monte::OccLocation const &occ_location,
              std::optional<std::vector<EventFilterGroup>> _event_filters,
              std::shared_ptr<engine_type> engine) override;

  /// \brief Run the KMC simulation
  void run(state_type &state, monte::OccLocation &occ_location,
           kmc_data_type &kmc_data, SelectedEvent &selected_event,
           std::optional<monte::SelectedEventDataCollector> &collector,
           run_manager_type &run_manager,
           std::shared_ptr<occ_events::OccSystem> event_system) override;

  /// \brief Return event selector type name
  std::string event_selector_type_str() const;

  /// \brief Constructs `event_selector` from the current `event_list` and
  /// `random_generator`; must be called after `update`
  void make_event_selector();

  // -- Validators --

  EventStateCalculator const &_prim_event_calculator(
      Index prim_event_index) const {
    if (prim_event_calculators.size() == 0) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData: "
          "prim_event_calculators.size() == 0");
    }
    if (prim_event_index >= prim_event_calculators.size()) {
      throw std::runtime_error(
          "AllowedKineticEventData::kra_coefficients: "
          "prim_event_index (=" +
          std::to_string(prim_event_index) +
          ") >= prim_event_calculators.size()");
    }
    return prim_event_calculators.at(prim_event_index);
  }

  AllowedEventList const &_event_list() const {
    if (!event_list) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData: Event list not set");
    }
    return *event_list;
  }

  Index _event_index(EventID const &id) const {
    AllowedEventList const &list = _event_list();
    auto it = list.allowed_event_map.find(id);
    if (it == list.allowed_event_map.events().end()) {
      throw std::runtime_error(
          "AllowedKineticEventData: Event not found in event list");
    }
    return std::distance(list.allowed_event_map.events().begin(), it);
  }

  AllowedEventCalculator &_event_calculator() const {
    if (!event_calculator) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData: Event calculator not set");
    }
    return *event_calculator;
  }

  event_selector_type const &_event_selector() const {
    if (!event_selector) {
      throw std::runtime_error(
          "Error in AllowedKineticEventData: Event selector not set");
    }
    return *event_selector;
  }

  // --- BaseMonteEventData interface ---

  // -- System data --

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients()
      const override {
    return _prim_event_calculator(0).formation_energy_coefficients();
  }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const override {
    return _prim_event_calculator(prim_event_index).freq_coefficients();
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const override {
    return _prim_event_calculator(prim_event_index).kra_coefficients();
  }

  // --- Event selection ---

  /// \brief Select an event, and optionally re-calculate event state for the
  ///     selected event
  void select_event(SelectedEvent &selected_event,
                    bool requires_event_state) override;

  /// Return number of events calculated with no barrier, by type
  std::map<std::string, Index> const &n_not_normal() const override {
    return _event_calculator().n_not_normal;
  }

  // -- Event list summary info --

  /// The size of the event list
  Index n_events() const override {
    return _event_list().allowed_event_map.n_assigned();
  }

  /// Return the current total event rate
  double total_rate() const override { return _event_selector().total_rate(); }

  // -- Event list iteration --

  /// Construct new internal iterator and return its index
  Index new_iterator(bool is_end) override {
    Index i = 0;
    while (m_it.find(i) != m_it.end()) {
      ++i;
    }
    auto const &x = _event_list().allowed_event_map.events();
    if (is_end) {
      m_it.emplace(i, x.end());
    } else {
      auto it = x.begin();
      while (it != x.end() && !it->is_assigned) {
        ++it;
      }
      m_it.emplace(i, it);
    }
    return i;
  }

  /// Copy internal iterator and return the new iterator index
  Index copy_iterator(Index i) override {
    if (m_it.find(i) == m_it.end()) {
      throw std::runtime_error(
          "AllowedKineticEventData::copy_iterator: Iterator not found");
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
          "AllowedKineticEventData::equal_iterator: Iterator not found");
    }
    return it_i->second == it_j->second;
  }

  /// Advance internal iterator by one event
  void advance_iterator(Index i) override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "AllowedKineticEventData::advance_iterator: Iterator not found");
    }
    auto &bare_it = it->second;
    auto end = _event_list().allowed_event_map.events().end();
    if (bare_it == end) {
      throw std::runtime_error(
          "AllowedKineticEventData::advance_iterator: Cannot advance past end "
          "of event list");
    }
    do {
      ++bare_it;
    } while (bare_it != end && !bare_it->is_assigned);
  }

  /// The event ID for the current state of the internal iterator
  EventID const &event_id(Index i) const override {
    auto it = m_it.find(i);
    if (it == m_it.end()) {
      throw std::runtime_error(
          "AllowedKineticEventData::event_id: Iterator not found");
    }
    return it->second->event_id;
  }

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the specified event. Reference is
  /// valid until the next call to this method.
  monte::OccEvent const &event_to_apply(EventID const &id) const override {
    return _event_calculator().set_event_data(id).event;
  }

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const override {
    auto it = _event_list().allowed_event_map.find(id);
    auto end = _event_list().allowed_event_map.events().end();
    if (it == end) {
      return 0.0;
    }
    auto begin = _event_list().allowed_event_map.events().begin();
    Index index = std::distance(begin, it);

    return _event_selector().get_rate(index);
  }

  /// Calculate event state data
  ///
  /// Notes:
  /// - Event state is only valid for event `id` until the next call to this
  ///   method
  EventState const &event_state(EventID const &id) const override {
    _event_calculator().calculate_rate(id);
    return _event_calculator().event_state;
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &impact(EventID const &id) const override {
    auto const &list = _event_list();
    return list.use_neighborlist_impact_table
               ? list.neighborlist_impact_table.value()(id)
               : list.relative_impact_table.value()(id);
  }

 private:
  /// \brief Holds internal iterators to allow generic interface for
  /// iterating over EventID
  std::map<Index, std::vector<AllowedEventData>::const_iterator> m_it;

  /// \brief Holds temporary calculated event state
  mutable EventState m_event_state;

  /// \brief Holds temporary calculated impact list
  std::vector<EventID> m_impact_list;
};

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

#endif