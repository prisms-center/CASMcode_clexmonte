#include "casm/clexmonte/monte_calculator/kinetic_events.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/events/io/stream/EventState_stream_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/configuration/Configuration.hh"
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

/// \brief Calculate the composition by orbit around an event and set
///     EventState::num_each_component_by_orbit
void LocalOrbitCompositionCalculator::calculate_num_each_component(
    EventState &state, EventData const &event_data,
    PrimEventData const &prim_event_data) {
  std::vector<Index> const &neighbor_index_to_linear_site_index =
      m_supercell_nlist->sites(event_data.unitcell_index);
  Eigen::VectorXi const &occupation = get_occupation(*m_state);

  // indices[orbit_index] = std::set<std::pair<int, int>>
  std::vector<std::set<std::pair<int, int>>> const &indices =
      m_local_orbits_neighbor_indices[prim_event_data.equivalent_index];

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

  state.num_each_component_by_orbit = &m_num_each_component_by_orbit;
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
    EventState &state, EventData const &event_data,
    PrimEventData const &prim_event_data) const {
  clexulator::ConfigDoFValues const *dof_values =
      m_formation_energy_clex->get();

  // Check if event is allowed based on the current occupation
  int i = 0;
  for (Index l : event_data.event.linear_site_index) {
    if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
      state.is_allowed = false;
      state.rate = 0.0;
      return;
    }
    ++i;
  }
  state.is_allowed = true;

  // calculate change in energy to final state
  //  state.dE_final = m_formation_energy_clex->occ_delta_value(
  //      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate change in energy to final state
  // - and save pointer to delta correlations
  state.formation_energy_delta_corr =
      &m_formation_energy_clex->correlations().occ_delta(
          event_data.event.linear_site_index, prim_event_data.occ_final);
  state.dE_final = m_formation_energy_clex->coefficients() *
                   (*state.formation_energy_delta_corr);

  // calculate KRA and attempt frequency
  // - add save pointer to local correlations
  state.local_corr = &m_event_clex->correlations().local(
      event_data.unitcell_index, prim_event_data.equivalent_index);
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

CompleteEventCalculator::CompleteEventCalculator(
    std::vector<PrimEventData> const &_prim_event_list,
    std::vector<EventStateCalculator> const &_prim_event_calculators,
    std::map<EventID, EventData> const &_event_list, Log &_event_log)
    : prim_event_list(_prim_event_list),
      prim_event_calculators(_prim_event_calculators),
      event_list(_event_list),
      event_log(_event_log),
      not_normal_count(0) {}

/// \brief Get CASM::monte::OccEvent corresponding to given event ID
double CompleteEventCalculator::calculate_rate(EventID const &id) {
  EventData const &event_data = event_list.at(id);
  PrimEventData const &prim_event_data =
      prim_event_list.at(id.prim_event_index);
  // Note: to keep all event state calculations, uncomment this:
  // EventState &event_state = event_data.event_state;
  prim_event_calculators.at(id.prim_event_index)
      .calculate_event_state(event_state, event_data, prim_event_data);

  // ---
  // can check event state and handle non-normal event states here
  // ---
  if (event_state.is_allowed && !event_state.is_normal) {
    event_log << "---" << std::endl;
    print(event_log.ostream(), event_state, event_data, prim_event_data);
    event_log << std::endl;
    ++not_normal_count;
  }

  return event_state.rate;
}

KineticEventData::KineticEventData(
    std::shared_ptr<system_type> _system,
    std::optional<std::vector<EventFilterGroup>> _event_filters)
    : transformation_matrix_to_super(Eigen::Matrix3l::Zero(3, 3)) {
  system = _system;
  if (!is_clex_data(*system, "formation_energy")) {
    throw std::runtime_error(
        "Error constructing KineticEventData: no 'formation_energy' clex.");
  }

  prim_event_list = clexmonte::make_prim_event_list(*system);
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
void KineticEventData::update(
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

    // TODO: rejection-clexmonte option does not require impact table
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
  event_selector = std::make_shared<KineticEventData::event_selector_type>(
      event_calculator,
      clexmonte::make_complete_event_id_list(n_unitcells, prim_event_list),
      event_list.impact_table,
      std::make_shared<lotto::RandomGenerator>(engine));
}

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM
