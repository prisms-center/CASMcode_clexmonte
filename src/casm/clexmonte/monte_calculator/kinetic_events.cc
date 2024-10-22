#include "casm/clexmonte/monte_calculator/kinetic_events.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

LocalOrbitCompositionCalculator::LocalOrbitCompositionCalculator(
    std::shared_ptr<system_type> _system, std::string _event_type_name,
    std::set<int> _orbits_to_calculate)
    : m_system(_system),
      m_event_type_name(_event_type_name),
      m_orbits_to_calculate(_orbits_to_calculate) {
  // Make m_occ_index_to_component_index_converter
  auto const &composition_calculator = get_composition_calculator(*m_system);
  m_occ_index_to_component_index_converter =
      composition::make_occ_index_to_component_index_converter(
          composition_calculator.components(),
          composition_calculator.allowed_occs());

  // Setup m_num_each_component_by_orbit and validate orbits_to_calculate
  auto cluster_info =
      get_local_basis_set_cluster_info(*m_system, m_event_type_name);
  int n_orbits = 0;
  if (cluster_info->orbits.size() > 0) {
    n_orbits = cluster_info->orbits[0].size();
  }

  for (int orbit_index : m_orbits_to_calculate) {
    if (orbit_index < 0 || orbit_index >= n_orbits) {
      std::stringstream msg;
      msg << "Error in LocalOrbitCompositionCalculator: "
          << "orbit_to_calculate=" << orbit_index << " out of range [0,"
          << n_orbits << ").";
      throw std::runtime_error(msg.str());
    }
  }

  m_num_each_component_by_orbit.resize(
      composition_calculator.components().size(), m_orbits_to_calculate.size());
  m_num_each_component_by_orbit.setZero();
}

/// \brief Reset pointer to state currently being calculated
void LocalOrbitCompositionCalculator::set(state_type const *state) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting LocalOrbitCompositionCalculator state: state is "
        "empty");
  }

  // set shell composition calculation data
  auto cluster_info =
      get_local_basis_set_cluster_info(*m_system, m_event_type_name);

  // Make m_local_orbits_neighbor_indices:
  m_local_orbits_neighbor_indices.clear();
  m_supercell_nlist = get_supercell_neighbor_list(*m_system, *m_state);
  auto const &convert = get_index_conversions(*m_system, *m_state);
  auto const &supercell_index_converter = convert.index_converter();
  for (Index equivalent_index = 0;
       equivalent_index < cluster_info->orbits.size(); ++equivalent_index) {
    std::vector<std::set<std::pair<int, int>>> _neighbor_indices_by_orbit;
    for (auto const &orbit : cluster_info->orbits[equivalent_index]) {
      std::set<std::pair<int, int>> _neighbor_indices;
      for (auto const &cluster : orbit) {
        for (auto const &site : cluster) {
          Index site_index = supercell_index_converter(site);
          _neighbor_indices.emplace(
              m_supercell_nlist->neighbor_index(site_index), site.sublattice());
        }
      }
      _neighbor_indices_by_orbit.emplace_back(std::move(_neighbor_indices));
    }
    m_local_orbits_neighbor_indices.emplace_back(
        std::move(_neighbor_indices_by_orbit));
  }
}

/// \brief Calculate the composition by orbit around an event
Eigen::MatrixXi const &
LocalOrbitCompositionCalculator::calculate_num_each_component(
    Eigen::VectorXi const &occupation, Index unitcell_index,
    Index equivalent_index) {
  std::vector<Index> const &neighbor_index_to_linear_site_index =
      m_supercell_nlist->sites(unitcell_index);

  // indices[orbit_index] = std::set<std::pair<int, int>>
  std::vector<std::set<std::pair<int, int>>> const &indices =
      m_local_orbits_neighbor_indices[equivalent_index];

  m_num_each_component_by_orbit.setZero();
  int col = 0;
  for (int orbit_index : m_orbits_to_calculate) {
    for (auto const &pair : indices[orbit_index]) {
      int neighbor_index = pair.first;
      int sublattice_index = pair.second;
      int site_index = neighbor_index_to_linear_site_index[neighbor_index];
      int occ_index = occupation(site_index);
      int component_index =
          m_occ_index_to_component_index_converter[sublattice_index][occ_index];
      m_num_each_component_by_orbit.col(col)(component_index) += 1;
    }
    ++col;
  }

  return m_num_each_component_by_orbit;
}

/// \brief Constructor
EventStateCalculator::EventStateCalculator(std::shared_ptr<system_type> _system,
                                           std::string _event_type_name)
    : m_system(_system), m_event_type_name(_event_type_name) {}

/// \brief Reset pointer to state currently being calculated
void EventStateCalculator::set(state_type const *state) {
  // supercell-specific
  m_state = state;
  if (m_state == nullptr) {
    throw std::runtime_error(
        "Error setting EventStateCalculator state: state is empty");
  }
  m_temperature = &m_state->conditions.scalar_values.at("temperature");
  m_formation_energy_clex = get_clex(*m_system, *m_state, "formation_energy");

  // set and validate event clex
  LocalMultiClexData event_local_multiclex_data =
      get_local_multiclex_data(*m_system, m_event_type_name);
  m_event_clex = get_local_multiclex(*m_system, *m_state, m_event_type_name);
  m_event_values.resize(m_event_clex->coefficients().size());
  std::map<std::string, Index> _glossary =
      event_local_multiclex_data.coefficients_glossary;

  auto _check_coeffs = [&](Index &coeff_index, std::string key) {
    if (!_glossary.count(key)) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: No " << key << " cluster expansion";
      throw std::runtime_error(ss.str());
    }
    coeff_index = _glossary.at(key);
    if (coeff_index < 0 || coeff_index >= m_event_clex->coefficients().size()) {
      std::stringstream ss;
      ss << "Error constructing " << m_event_type_name
         << " EventStateCalculator: " << key << " index out of range";
      throw std::runtime_error(ss.str());
    }
  };
  _check_coeffs(m_kra_index, "kra");
  _check_coeffs(m_freq_index, "freq");
}

/// \brief Pointer to current state
state_type const *EventStateCalculator::state() const { return m_state; }

/// \brief Current state's reciprocal temperature
double EventStateCalculator::beta() const {
  return 1.0 / (CASM::KB * *this->m_temperature);
}

/// \brief Calculate the state of an event
void EventStateCalculator::calculate_event_state(
    EventState &state, Index unitcell_index,
    std::vector<Index> const &linear_site_index,
    PrimEventData const &prim_event_data) const {
  clexulator::ConfigDoFValues const *dof_values =
      m_formation_energy_clex->get();

  // Check if event is allowed based on the current occupation
  state.is_allowed =
      event_is_allowed(linear_site_index, *dof_values, prim_event_data);
  if (!state.is_allowed) {
    state.rate = 0.0;
    return;
  }

  // calculate change in energy to final state
  //  state.dE_final = m_formation_energy_clex->occ_delta_value(
  //      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate change in energy to final state
  // - and save pointer to delta correlations
  state.formation_energy_delta_corr =
      &m_formation_energy_clex->correlations().occ_delta(
          linear_site_index, prim_event_data.occ_final);
  state.dE_final = m_formation_energy_clex->coefficients() *
                   (*state.formation_energy_delta_corr);

  // calculate KRA and attempt frequency
  // - add save pointer to local correlations
  state.local_corr = &m_event_clex->correlations().local(
      unitcell_index, prim_event_data.equivalent_index);
  for (int i = 0; i < m_event_clex->coefficients().size(); ++i) {
    m_event_values(i) = m_event_clex->coefficients()[i] * (*state.local_corr);
  }
  state.Ekra = m_event_values[m_kra_index];
  state.freq = m_event_values[m_freq_index];

  // calculate energy in activated state, check if "normal", calculate rate
  state.dE_activated = state.dE_final * 0.5 + state.Ekra;
  state.is_normal =
      (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
  if (state.dE_activated < state.dE_final) state.dE_activated = state.dE_final;
  if (state.dE_activated < 0.0) state.dE_activated = 0.0;
  state.rate = state.freq * exp(-this->beta() * state.dE_activated);
}

namespace {

void print_no_barrier_warning(Log &event_log, EventState const &event_state,
                              EventData const &event_data,
                              PrimEventData const &prim_event_data) {
  event_log << "## WARNING: EVENT WITH NO BARRIER ###################\n"
               "#                                                   #\n"
               "# Events with no barrier are treated as having a    #\n"
               "# rate equal to the attempt frequency.              #\n"
               "#                                                   #\n"
               "# This warning is only printed once per event type. #\n"
               "#                                                   #\n"
               "# Event info:                                       #\n"
            << std::endl;
  print(event_log.ostream(), event_state, event_data, prim_event_data);
  event_log << "#                                                   #\n"
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

CompleteKineticEventData::CompleteKineticEventData(
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
void CompleteKineticEventData::update(
    state_type const &state, monte::OccLocation const &occ_location,
    std::optional<std::vector<EventFilterGroup>> _event_filters,
    std::shared_ptr<engine_type> engine) {
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

/// \brief Update for given state, conditions, occupants, event filters
void CompleteKineticEventData::select_event(SelectedEvent &selected_event,
                                            bool requires_event_state) {
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

AllowedKineticEventData::AllowedKineticEventData(
    std::shared_ptr<system_type> _system, bool _allow_events_with_no_barrier,
    bool _use_map_index, bool _use_neighborlist_impact_table,
    kinetic_event_selector_type _event_selector_type)
    : allow_events_with_no_barrier(_allow_events_with_no_barrier),
      use_map_index(_use_map_index),
      use_neighborlist_impact_table(_use_neighborlist_impact_table),
      event_selector_type(_event_selector_type) {
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

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});
}

/// \brief Update for given state, conditions, occupants, event filters
///
/// Notes:
/// - This constructs the complete event list and impact table, and constructs
///   the event selector, which calculates all event rates.
/// - If there are no event filters and the supercell remains unchanged from the
///   previous update, then the event list and impact table are not
///   re-constructed, but the event rates are still re-calculated.
void AllowedKineticEventData::update(state_type const &state,
                                     monte::OccLocation const &occ_location,
                                     std::shared_ptr<engine_type> engine) {
  random_generator = std::make_shared<lotto::RandomGenerator>(engine);

  // These are constructed/re-constructed so cluster expansions point
  // at the current state
  prim_event_calculators.clear();
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(system,
                                        prim_event_data.event_type_name);
    prim_event_calculators.back().set(&state);
  }

  // Construct AllowedEventList
  event_list = std::make_shared<clexmonte::AllowedEventList>(
      prim_event_list, prim_impact_info_list, get_dof_values(state),
      occ_location, get_prim_neighbor_list(*system),
      get_supercell_neighbor_list(*system, state), use_map_index,
      use_neighborlist_impact_table);

  if (event_list->allowed_event_map.n_assigned() == 0) {
    throw std::runtime_error(
        "Error constructing event list: "
        "no allowed events.");
  }

  // Construct AllowedEventCalculator
  event_calculator = std::make_shared<AllowedEventCalculator>(
      prim_event_list, prim_event_calculators, *event_list);

  // Make event selector
  // - This calculates all rates at construction
  make_event_selector();
}

/// \brief Constructs `event_selector` from the current `event_list` and
/// `random_generator`; must be called after `update`
void AllowedKineticEventData::make_event_selector() {
  // Make event selector
  // - This calculates all rates at construction

  if (this->event_selector_type == kinetic_event_selector_type::sum_tree) {
    std::cout << "Constructing \"sum_tree\" event selector" << std::endl;
    sum_tree_event_selector = std::make_shared<sum_tree_event_selector_type>(
        event_calculator,
        this->event_list->allowed_event_map.event_index_list(),
        GetImpactFromAllowedEventList(this->event_list),
        this->random_generator);
  } else if (this->event_selector_type ==
             kinetic_event_selector_type::vector_sum_tree) {
    std::cout << "Constructing \"vector_sum_tree\" event selector" << std::endl;
    vector_sum_tree_event_selector =
        std::make_shared<vector_sum_tree_event_selector_type>(
            event_calculator,
            this->event_list->allowed_event_map.events().size(),
            GetImpactFromAllowedEventList(this->event_list),
            this->random_generator);
  } else {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: "
        "invalid event_selector_type.");
  }
}

/// \brief Update for given state, conditions, occupants, event filters
void AllowedKineticEventData::select_event(SelectedEvent &selected_event,
                                           bool requires_event_state) {
  if (this->event_list->allowed_event_map.has_new_events()) {
    this->make_event_selector();
    this->event_list->allowed_event_map.clear_has_new_events();
  }

  // This function:
  // - Updates rates of events impacted by the *last* selected event (if there
  //   was a previous selection)
  // - Updates the total rate
  // - Chooses an event and time increment (does not apply event)
  // - Sets a list of impacted events by the chosen event
  Index selected_event_index;

  if (this->event_selector_type == kinetic_event_selector_type::sum_tree) {
    // sum_tree_event_selector
    std::tie(selected_event_index, selected_event.time_increment) =
        sum_tree_event_selector->select_event();
    selected_event.total_rate = sum_tree_event_selector->total_rate();

  } else if (this->event_selector_type ==
             kinetic_event_selector_type::vector_sum_tree) {
    // vector_sum_tree_event_selector
    std::tie(selected_event_index, selected_event.time_increment) =
        vector_sum_tree_event_selector->select_event();
    selected_event.total_rate = vector_sum_tree_event_selector->total_rate();

  } else {
    throw std::runtime_error(
        "Error constructing AllowedKineticEventData: "
        "invalid event_selector_type.");
  }

  EventID const &event_id =
      this->event_list->allowed_event_map.event_id(selected_event_index);
  EventData const &event_data =
      event_calculator->set_event_data(selected_event_index);

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

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM
