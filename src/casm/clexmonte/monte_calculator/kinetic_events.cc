#include "casm/clexmonte/monte_calculator/kinetic_events.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/State.hh"

// this must be included after the classes are defined for its application here
#include "casm/clexmonte/methods/kinetic_monte_carlo.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

namespace {

void print_no_barrier_warning(Log &event_log, EventState const &event_state,
                              EventData const &event_data,
                              PrimEventData const &prim_event_data) {
  event_log << "\n"
               "## WARNING: EVENT WITH NO BARRIER ###################\n"
               "#                                                   #\n"
               "# Events with no barrier are treated as having a    #\n"
               "# rate equal to the attempt frequency.              #\n"
               "#                                                   #\n"
               "# This warning is only printed once per event type. #\n"
               "#                                                   #\n"
               "# Event info:                                       #\n"
            << std::endl;
  print(event_log, event_state, event_data, prim_event_data);
  event_log << "\n"
               "#                                                   #\n"
               "#####################################################\n"
            << std::endl;
}

}  // namespace

// -- CompleteKineticEventData --

CompleteEventCalculator::CompleteEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    std::map<EventID, EventData> const &_event_list, Log &_event_log)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_log(_event_log) {}

/// \brief Update `event_state` for event `id` in the current state and
/// return the event rate
double CompleteEventCalculator::calculate_rate(EventID const &id) {
  EventData const &event_data = event_list.at(id);
  PrimEventData const &prim_event_data =
      prim_event_list.at(id.prim_event_index);
  // Note: to keep all event state calculations, uncomment this:
  // EventState &event_state = event_data.event_state;
  prim_event_calculators.at(id.prim_event_index)
      .calculate_event_state(event_state, event_data.unitcell_index,
                             event_data.event.linear_site_index,
                             prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (event_state.is_allowed && !event_state.is_normal) {
    Index &n = n_not_normal[prim_event_data.event_type_name];

    if (n == 0) {
      print_no_barrier_warning(event_log, event_state, event_data,
                               prim_event_data);
    }
    n += 1;
  }

  return event_state.rate;
}

template <bool DebugMode>
CompleteKineticEventData<DebugMode>::CompleteKineticEventData(
    std::shared_ptr<system_type> _system,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    bool _allow_events_with_no_barrier)
    : allow_events_with_no_barrier(_allow_events_with_no_barrier),
      transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)) {
  system = _system;
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing CompleteKineticEventData: no 'formation_energy' "
        "clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
  if (prim_event_list.empty()) {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: "
        "prim event list is empty.");
  }

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});

  if (_event_filters.has_value()) {
    event_filters = _event_filters.value();
  }
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// Notes:
/// - This constructs the complete event list and impact table, and constructs
///   the event selector, which calculates all event rates.
/// - If there are no event filters and the supercell remains unchanged from the
///   previous update, then the event list and impact table are not
///   re-constructed, but the event rates are still re-calculated.
template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::update(
    std::shared_ptr<StateData> _state_data,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    std::shared_ptr<engine_type> engine) {
  // Current state info
  state_data = _state_data;
  state_type const &state = *state_data->state;
  monte::OccLocation const &occ_location = *state_data->occ_location;

  // if same supercell && no event filters
  // -> just re-set state & avoid re-constructing event list
  if (this->transformation_matrix_to_super ==
          get_transformation_matrix_to_super(state) &&
      !_event_filters.has_value()) {
    for (auto &event_state_calculator : prim_event_calculators) {
      event_state_calculator.set(&state);
    }
    event_calculator->n_not_normal.clear();
  } else {
    if (_event_filters.has_value()) {
      event_filters = _event_filters.value();
    }

    // These are constructed/re-constructed so cluster expansions point
    // at the current state
    prim_event_calculators.clear();
    for (auto const &prim_event_data : prim_event_list) {
      prim_event_calculators.emplace_back(system,
                                          prim_event_data.event_type_name);
      prim_event_calculators.back().set(&state);

      // Set a custom event state calculation function if it exists:
      auto it = custom_event_state_calculation_f.find(
          prim_event_data.event_type_name);
      if (it != custom_event_state_calculation_f.end()) {
        prim_event_calculators.back().set_custom_event_state_calculation(
            it->second);
      }
    }

    // Construct CompleteEventList
    event_list = clexmonte::make_complete_event_list(
        prim_event_list, prim_impact_info_list, occ_location, event_filters);

    // Construct CompleteEventCalculator
    event_calculator = std::make_shared<CompleteEventCalculator>(
        prim_event_list, prim_event_calculators, event_list.events);

    transformation_matrix_to_super = get_transformation_matrix_to_super(state);
  }

  Index n_unitcells = transformation_matrix_to_super.determinant();

  // Make event selector
  // - This calculates all rates at construction
  event_selector =
      std::make_shared<CompleteKineticEventData::event_selector_type>(
          event_calculator,
          clexmonte::make_complete_event_id_list(n_unitcells, prim_event_list),
          event_list.impact_table,
          std::make_shared<lotto::RandomGenerator>(engine));
}

template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::run(
    state_type &state, monte::OccLocation &occ_location,
    kmc_data_type &kmc_data, SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Function to set selected event
  bool requires_event_state =
      collector.has_value() && collector->requires_event_state;
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->select_event(selected_event, requires_event_state);
  };

  auto set_impacted_events_f = [=](SelectedEvent &selected_event) {
    // Set impacted events
    this->event_selector->set_impacted_events(selected_event.event_id);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, kmc_data, selected_event, set_selected_event_f,
      set_impacted_events_f, collector, run_manager, event_system);
}

/// \brief Update for given state, conditions, occupants, event filters
template <bool DebugMode>
void CompleteKineticEventData<DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  // This function:
  // - Updates rates of events impacted by the *last* selected event (if there
  //   was a previous selection)
  // - Updates the total rate
  // - Chooses an event and time increment (does not apply event)
  // - Sets a list of impacted events by the chosen event
  std::tie(selected_event.event_id, selected_event.time_increment) =
      event_selector->select_event();
  selected_event.total_rate = event_selector->total_rate();
  EventID const &event_id = selected_event.event_id;
  EventData const &event_data = event_list.events.at(event_id);
  PrimEventData const &prim_event_data =
      prim_event_list[event_id.prim_event_index];
  selected_event.event_data = &event_data;
  selected_event.prim_event_data = &prim_event_data;

  if (!allow_events_with_no_barrier && event_calculator->n_not_normal.size()) {
    throw std::runtime_error(
        "Error: Encountered event with no barrier, which is not allowed.");
  }

  if (requires_event_state) {
    prim_event_calculators.at(event_id.prim_event_index)
        .calculate_event_state(m_event_state, event_data.unitcell_index,
                               event_data.event.linear_site_index,
                               prim_event_data);
    selected_event.event_state = &m_event_state;
  }
}

// -- AllowedKineticEventData --

AllowedEventCalculator::AllowedEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    AllowedEventList &_event_list, Log &_event_log)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_log(_event_log) {}

/// \brief Update `event_state` for event `event_index` in the current state
/// and return the event rate; if the event is no longer allowed, free the
/// event.
double AllowedEventCalculator::calculate_rate(Index event_index) {
  AllowedEventData const &allowed_event_data =
      event_list.allowed_event_map.events()[event_index];
  // EventID original_event_id = allowed_event_data.event_id;
  if (!allowed_event_data.is_assigned) {
    event_state.is_allowed = false;
    event_state.rate = 0.0;
  } else {
    this->calculate_rate(allowed_event_data.event_id);

    // free event from AllowedEventList if not allowed
    if (!event_state.is_allowed) {
      event_list.allowed_event_map.free(allowed_event_data.event_id);
    }
  }

  return event_state.rate;
}

/// \brief Update `event_state` for any event `event_id` in the current state
/// and return the event rate
double AllowedEventCalculator::calculate_rate(EventID const &event_id) {
  Index prim_event_index = event_id.prim_event_index;
  PrimEventData const &prim_event_data =
      this->prim_event_list[prim_event_index];
  Index unitcell_index = event_id.unitcell_index;

  // set linear_site_index
  set_event_linear_site_index(linear_site_index, unitcell_index,
                              event_list.neighbor_index[prim_event_index],
                              *event_list.supercell_nlist);

  // calculate event state
  prim_event_calculators.at(prim_event_index)
      .calculate_event_state(event_state, unitcell_index, linear_site_index,
                             prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (event_state.is_allowed && !event_state.is_normal) {
    Index &n = n_not_normal[prim_event_data.event_type_name];

    if (n == 0) {
      set_event_data(event_id);
      print_no_barrier_warning(event_log, event_state, event_data,
                               prim_event_data);
    }
    n += 1;
  }

  return event_state.rate;
}

/// \brief Set `event_data` for event `event_index`, returning a reference
/// which is valid until the next call to this method
EventData const &AllowedEventCalculator::set_event_data(Index event_index) {
  return set_event_data(event_list.allowed_event_map.event_id(event_index));
}

/// \brief Set `event_data` for any event `event_id`, returning a reference
/// which is valid until the next call to this method
EventData const &AllowedEventCalculator::set_event_data(
    EventID const &event_id) {
  Index prim_event_index = event_id.prim_event_index;
  PrimEventData const &prim_event_data =
      this->prim_event_list[prim_event_index];
  Index unitcell_index = event_id.unitcell_index;

  // set this->event_data.unitcell_index
  this->event_data.unitcell_index = unitcell_index;

  // set this->event_data.event
  set_event(this->event_data.event, prim_event_data, unitcell_index,
            event_list.occ_location,
            event_list.neighbor_index[prim_event_index],
            *event_list.supercell_nlist);

  return this->event_data;
}

template <typename EventSelectorType, bool DebugMode>
AllowedKineticEventData<EventSelectorType, DebugMode>::AllowedKineticEventData(
    std::shared_ptr<system_type> _system, bool _allow_events_with_no_barrier,
    bool _use_map_index, bool _use_neighborlist_impact_table,
    bool _assign_allowed_events_only)
    : allow_events_with_no_barrier(_allow_events_with_no_barrier),
      use_map_index(_use_map_index),
      use_neighborlist_impact_table(_use_neighborlist_impact_table),
      assign_allowed_events_only(_assign_allowed_events_only) {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Construct AllowedKineticEventData");
    log << "Event data and selection:" << std::endl;
    log << "- impact_table_type="
        << (use_neighborlist_impact_table ? std::string("\"neighborlist\"")
                                          : std::string("\"relative\""))
        << std::endl;
    log << "- event_selector_type=\"" << this->event_selector_type_str() << "\""
        << std::endl;
    log << "- assign_allowed_events_only=" << std::boolalpha
        << assign_allowed_events_only << std::endl;
    log << std::endl;
    log.end_section();
  }

  system = _system;
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: no 'formation_energy' "
        "clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
  if (prim_event_list.empty()) {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: "
        "prim event list is empty.");
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Prim event list");
    log.indent() << qto_json(prim_event_list) << std::endl << std::endl;
  }

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// Notes:
/// - This constructs the complete event list and impact table, and constructs
///   the event selector, which calculates all event rates.
/// - Event filters are ignored (with a warning). This is a TODO feature.
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::update(
    std::shared_ptr<StateData> _state_data,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    std::shared_ptr<engine_type> engine) {
  random_generator = std::make_shared<lotto::RandomGenerator>(engine);
  state_data = _state_data;

  // Warning if event_filters:
  if (_event_filters.has_value()) {
    std::cerr << "#############################################" << std::endl;
    std::cerr << "Warning: Event filters are being ignored. Use" << std::endl;
    std::cerr << "the \"high_memory\" event data type to apply " << std::endl;
    std::cerr << "event filters.                               " << std::endl;
    std::cerr << "#############################################" << std::endl;
  }

  // Current state info
  state_type const &state = *state_data->state;
  monte::OccLocation const &occ_location = *state_data->occ_location;

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Monte Carlo State");
    log.indent() << qto_json(state) << std::endl << std::endl;
  }

  // These are constructed/re-constructed so cluster expansions point
  // at the current state
  prim_event_calculators.clear();
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(system,
                                        prim_event_data.event_type_name);
    prim_event_calculators.back().set(&state);

    // Set a custom event state calculation function if it exists:
    auto it =
        custom_event_state_calculation_f.find(prim_event_data.event_type_name);
    if (it != custom_event_state_calculation_f.end()) {
      prim_event_calculators.back().set_custom_event_state_calculation(
          it->second);
    }
  }

  // Construct AllowedEventList
  event_list = std::make_shared<clexmonte::AllowedEventList>(
      prim_event_list, prim_impact_info_list, get_dof_values(state),
      occ_location, get_prim_neighbor_list(*system),
      get_supercell_neighbor_list(*system, state), use_map_index,
      use_neighborlist_impact_table, assign_allowed_events_only);

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Event list summary");
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log << std::endl;
    log.end_section();
  }

  // Construct AllowedEventCalculator
  event_calculator = std::make_shared<AllowedEventCalculator>(
      prim_event_list, prim_event_calculators, *event_list);

  // Make event selector
  // - This calculates all rates at construction
  this->make_event_selector();
  this->event_list->allowed_event_map.clear_has_new_events();
}

template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::run(
    state_type &state, monte::OccLocation &occ_location,
    kmc_data_type &kmc_data, SelectedEvent &selected_event,
    std::optional<monte::SelectedEventDataCollector> &collector,
    run_manager_type &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Function to set selected event
  bool requires_event_state =
      collector.has_value() && collector->requires_event_state;
  auto set_selected_event_f = [=](SelectedEvent &selected_event) {
    this->select_event(selected_event, requires_event_state);
  };

  auto set_impacted_events_f = [=](SelectedEvent &selected_event) {
    // Set impacted events
    this->event_selector->set_impacted_events(selected_event.event_index);
  };

  // Run Kinetic Monte Carlo at a single condition
  kinetic_monte_carlo_v2<DebugMode>(
      state, occ_location, kmc_data, selected_event, set_selected_event_f,
      set_impacted_events_f, collector, run_manager, event_system);
}

// -- EventSelectorType specializations --
namespace {
template <typename EventSelectorType>
std::string event_selector_type_str_impl();

template <>
std::string event_selector_type_str_impl<sum_tree_event_selector_type>() {
  return "sum_tree";
}

template <>
std::string
event_selector_type_str_impl<vector_sum_tree_event_selector_type>() {
  return "vector_sum_tree";
}

template <>
std::string event_selector_type_str_impl<direct_sum_event_selector_type>() {
  return "direct_sum";
}

template <typename EventSelectorType>
std::shared_ptr<EventSelectorType> make_event_selector_impl(
    std::shared_ptr<AllowedEventCalculator> event_calculator,
    std::shared_ptr<AllowedEventList> event_list,
    std::shared_ptr<lotto::RandomGenerator> random_generator);

template <>
std::shared_ptr<sum_tree_event_selector_type> make_event_selector_impl(
    std::shared_ptr<AllowedEventCalculator> event_calculator,
    std::shared_ptr<AllowedEventList> event_list,
    std::shared_ptr<lotto::RandomGenerator> random_generator) {
  return std::make_shared<sum_tree_event_selector_type>(
      event_calculator, event_list->allowed_event_map.event_index_list(),
      GetImpactFromAllowedEventList(event_list), random_generator);
}

template <>
std::shared_ptr<vector_sum_tree_event_selector_type> make_event_selector_impl(
    std::shared_ptr<AllowedEventCalculator> event_calculator,
    std::shared_ptr<AllowedEventList> event_list,
    std::shared_ptr<lotto::RandomGenerator> random_generator) {
  return std::make_shared<vector_sum_tree_event_selector_type>(
      event_calculator, event_list->allowed_event_map.events().size(),
      GetImpactFromAllowedEventList(event_list), random_generator);
}

template <>
std::shared_ptr<direct_sum_event_selector_type> make_event_selector_impl(
    std::shared_ptr<AllowedEventCalculator> event_calculator,
    std::shared_ptr<AllowedEventList> event_list,
    std::shared_ptr<lotto::RandomGenerator> random_generator) {
  return std::make_shared<direct_sum_event_selector_type>(
      event_calculator, event_list->allowed_event_map.events().size(),
      GetImpactFromAllowedEventList(event_list), random_generator);
}

}  // namespace
// -- end EventSelectorType specializations --

template <typename EventSelectorType, bool DebugMode>
std::string AllowedKineticEventData<
    EventSelectorType, DebugMode>::event_selector_type_str() const {
  return event_selector_type_str_impl<EventSelectorType>();
}

/// \brief Constructs `event_selector`; must be called after `update`
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType,
                             DebugMode>::make_event_selector() {
  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Make event selector");
    log.indent() << "- event_selector_type=\""
                 << this->event_selector_type_str() << "\"" << std::endl;
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log.indent() << "- Constructing event selector..." << std::endl;
  }

  // Make event selector
  // - This calculates all rates at construction
  event_selector = make_event_selector_impl<EventSelectorType>(
      event_calculator, event_list, random_generator);

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.indent() << "- Constructing event selector... DONE" << std::endl;
    log.indent() << "- total_rate=" << this->event_selector->total_rate()
                 << std::endl;
    log << std::endl;
    log.end_section();
  }
}

/// \brief Update for given state, conditions, occupants, event filters
template <typename EventSelectorType, bool DebugMode>
void AllowedKineticEventData<EventSelectorType, DebugMode>::select_event(
    SelectedEvent &selected_event, bool requires_event_state) {
  // If updating the event list with impacted events after the previous step
  // caused the event list to increase in size, then it needs to be
  // re-constructed.
  if (this->event_list->allowed_event_map.has_new_events()) {
    if constexpr (DebugMode) {
      Log &log = CASM::log();
      log.custom("Select event requires re-constructing event selector");
      log << std::endl;
      CASM::log().increase_indent();
    }

    this->make_event_selector();
    this->event_list->allowed_event_map.clear_has_new_events();

    if constexpr (DebugMode) {
      CASM::log().decrease_indent();
    }
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Event list summary");
    log.indent() << "- Event list container size: "
                 << event_list->allowed_event_map.n_total() << std::endl;
    log.indent() << "- Number of events: "
                 << event_list->allowed_event_map.n_assigned() << std::endl;
    log << std::endl;
    log.end_section();
  }

  // The function `only_select_event` does the following:
  // - Updates rates of events impacted by the *last* selected event (if there
  //   was a previous selection)
  // - Updates the total rate
  // - Chooses an event and time increment
  //
  // It does not apply the event or set the impacted events.
  Index selected_event_index;
  std::tie(selected_event_index, selected_event.time_increment) =
      event_selector->only_select_event();
  selected_event.total_rate = event_selector->total_rate();

  EventID const &event_id =
      this->event_list->allowed_event_map.event_id(selected_event_index);
  EventData const &event_data =
      event_calculator->set_event_data(selected_event_index);
  PrimEventData const &prim_event_data =
      prim_event_list[event_id.prim_event_index];

  selected_event.event_id = event_id;
  selected_event.event_index = selected_event_index;
  selected_event.event_data = &event_data;
  selected_event.prim_event_data = &prim_event_data;

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.custom("Selected event");

    prim_event_calculators.at(event_id.prim_event_index)
        .calculate_event_state(m_event_state, event_data.unitcell_index,
                               event_data.event.linear_site_index,
                               prim_event_data);
    selected_event.event_state = &m_event_state;

    // get ijk
    auto const &unitcell_index_converter =
        state_data->occ_location->convert().unitcell_index_converter();
    auto ijk = unitcell_index_converter(event_id.unitcell_index);
    jsonParser ijk_json;
    to_json(ijk, ijk_json, jsonParser::as_array());

    // format output
    {
      jsonParser tjson;
      to_json(selected_event, tjson, *get_event_system(*system));

      // Add event position
      jsonParser ajson;
      ajson["unitcell_index"] = event_data.unitcell_index;
      ajson["unitcell_ijk"] = ijk_json;

      ajson["event_sites_relative"] = jsonParser::array();
      for (auto const &site : prim_event_data.sites) {
        ajson["event_sites_relative"].push_back(qto_json(site));
      }

      ajson["event_sites_absolute"] = jsonParser::array();
      for (auto const &site : prim_event_data.sites) {
        ajson["event_sites_absolute"].push_back(qto_json(site + ijk));
      }
      tjson["event_position"] = ajson;

      // Add Monte Carlo state
      tjson["state"] = *state_data->state;

      log << tjson << std::endl << std::endl;
    }
  }

  if (!allow_events_with_no_barrier && event_calculator->n_not_normal.size()) {
    throw std::runtime_error(
        "Error: Encountered event with no barrier, which is not allowed.");
  }

  if (requires_event_state) {
    prim_event_calculators.at(event_id.prim_event_index)
        .calculate_event_state(m_event_state, event_data.unitcell_index,
                               event_data.event.linear_site_index,
                               prim_event_data);
    selected_event.event_state = &m_event_state;
  }

  if constexpr (DebugMode) {
    Log &log = CASM::log();
    log.end_section();
  }
}

// Explicit instantiation:

// DebugMode=false
template class CompleteKineticEventData<false>;
template class AllowedKineticEventData<vector_sum_tree_event_selector_type,
                                       false>;
template class AllowedKineticEventData<sum_tree_event_selector_type, false>;
template class AllowedKineticEventData<direct_sum_event_selector_type, false>;

// DebugMode=true
template class CompleteKineticEventData<true>;
template class AllowedKineticEventData<vector_sum_tree_event_selector_type,
                                       true>;
template class AllowedKineticEventData<sum_tree_event_selector_type, true>;
template class AllowedKineticEventData<direct_sum_event_selector_type, true>;

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM
