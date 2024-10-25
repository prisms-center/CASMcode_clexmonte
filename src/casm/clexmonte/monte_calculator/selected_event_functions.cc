#include "casm/clexmonte/monte_calculator/selected_event_functions.hh"

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {
namespace monte_calculator {

namespace {
Index get_prim_event_index(
    std::shared_ptr<SelectedEvent> const &selected_event) {
  return selected_event->prim_event_data->prim_event_index;
}

Index get_equivalent_index(
    std::shared_ptr<SelectedEvent> const &selected_event) {
  return selected_event->prim_event_data->equivalent_index;
}

Index get_unitcell_index(std::shared_ptr<SelectedEvent> const &selected_event) {
  return selected_event->event_id.unitcell_index;
}

EventState const &get_event_state(
    std::shared_ptr<SelectedEvent> const &selected_event) {
  return *selected_event->event_state;
}

SelectedEventInfo make_selected_event_info(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return SelectedEventInfo(get_prim_event_list(calculation));
}

}  // namespace

/// \brief Make selected event count collecting function
/// ("selected_event.by_type")
monte::DiscreteVectorIntHistogramFunction make_selected_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  SelectedEventInfo info = make_selected_event_info(calculation);
  info.make_indices_by_type();
  auto prim_event_index_to_index = info.prim_event_index_to_index;
  auto selected_event = calculation->selected_event();

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_type",
      /* description */
      "Selected event count by event type",
      /* shape */
      {},
      /* component_names */
      std::nullopt,
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event]() {
        int value =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return to_VectorXl(value);
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = info.value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.by_equivalent_index")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  SelectedEventInfo info = make_selected_event_info(calculation);
  info.make_indices_by_equivalent_index();
  auto prim_event_index_to_index = info.prim_event_index_to_index;
  auto selected_event = calculation->selected_event();

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_equivalent_index",
      /* description */
      "Selected event count by event type and equivalent index",
      /* shape */
      {},
      /* component_names */
      std::nullopt,
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event]() {
        int value =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return to_VectorXl(value);
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = info.value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.by_equivalent_index_and_direction")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  SelectedEventInfo info = make_selected_event_info(calculation);
  info.make_indices_by_equivalent_index_and_direction();
  auto prim_event_index_to_index = info.prim_event_index_to_index;
  auto selected_event = calculation->selected_event();

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_equivalent_index_and_direction",
      /* description */
      "Selected event count by event type, equivalent index, and direction",
      /* shape */
      {},
      /* component_names */
      std::nullopt,
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event]() {
        int value =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return to_VectorXl(value);
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = info.value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.<event_type>.by_equivalent_index")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_selected_event_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // get names in alphabetical order
  std::set<std::string> keys;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    keys.insert(x.event_type_name);
  }

  std::vector<monte::DiscreteVectorIntHistogramFunction> f_list;

  for (auto const &event_type_name : keys) {
    SelectedEventInfo info = make_selected_event_info(calculation);
    info.make_indices_by_equivalent_index_per_event_type(event_type_name);
    auto prim_event_index_to_index = info.prim_event_index_to_index;
    auto prim_event_index_to_has_value = info.prim_event_index_to_has_value;
    auto selected_event = calculation->selected_event();

    // Create the function
    monte::DiscreteVectorIntHistogramFunction f(
        /* name */
        "selected_event." + event_type_name + ".by_equivalent_index",
        /* description */
        "Selected event count by equivalent index for event_type " +
            event_type_name,
        /* shape */
        {},
        /* component_names */
        std::nullopt,
        /* requires_event_state */
        false,
        /* function */
        [prim_event_index_to_index, selected_event]() {
          int value = (*prim_event_index_to_index)[get_prim_event_index(
              selected_event)];
          return to_VectorXl(value);
        },
        /* has_value_function */
        [prim_event_index_to_has_value, selected_event]() -> bool {
          return (*prim_event_index_to_has_value)[get_prim_event_index(
              selected_event)];
        },
        /* max_size */
        10000);
    f.value_labels = info.value_labels;
    f_list.push_back(f);
  }

  return f_list;
}

//  auto const &system = get_system(m_calculation);
//  auto const &composition_calculator = get_composition_calculator(system);
//  auto const &orbits =
//      get_local_basis_set_cluster_info(system, m_data.local_basis_set_name)
//          ->orbits;
//  std::shared_ptr<PossibleLocalOrbitCompositions const> compositions =
//      std::make_shared<PossibleLocalOrbitCompositions const>(
//          composition_calculator, orbits, m_data->max_size);

LocalOrbitCompositionCollector::LocalOrbitCompositionCollector(
    std::shared_ptr<MonteCalculator> calculation,
    std::shared_ptr<LocalOrbitCompositionCalculatorData const> data)
    : m_calculation(calculation), m_data(data) {}

/// \brief Return the shape of the output values
std::vector<Index> LocalOrbitCompositionCollector::shape() const {
  std::vector<Index> shape;
  auto const &composition_calculator =
      get_composition_calculator(get_system(m_calculation));
  shape.push_back(composition_calculator.components().size());
  if (m_data->combine_orbits) {
    shape.push_back(1);
  } else {
    shape.push_back(m_data->orbits_to_calculate.size());
  }
  return shape;
}

/// \brief Return names for value components (col-major unrolling)
std::vector<std::string> LocalOrbitCompositionCollector::component_names()
    const {
  std::vector<std::string> component_names;
  auto const &composition_calculator =
      get_composition_calculator(get_system(m_calculation));
  if (m_data->combine_orbits) {
    for (auto const &x : composition_calculator.components()) {
      component_names.push_back(x);
    }
  } else {
    for (int i_orbit : m_data->orbits_to_calculate) {
      for (auto const &x : composition_calculator.components()) {
        std::stringstream ss;
        ss << x << "," << i_orbit;
        component_names.push_back(ss.str());
      }
    }
  }
  return component_names;
}

/// \brief Evaluate the local orbit composition for a particular event
Eigen::MatrixXi const &LocalOrbitCompositionCollector::value(
    Index unitcell_index, Index equivalent_index) {
  if (!m_local_orbit_composition_calculator) {
    this->reset();
  }
  return m_local_orbit_composition_calculator->value(unitcell_index,
                                                     equivalent_index);
}

/// \brief Reset the LocalOrbitCompositionCalculator using the current state
/// being calculated by the MonteCalculator
void LocalOrbitCompositionCollector::reset() {
  auto &system = *m_calculation->system();
  if (!m_calculation->state_data()) {
    throw std::runtime_error(
        "Error in LocalOrbitCompositionCollector: state_data is "
        "not set in MonteCalculator");
  }
  auto const &state_data = *m_calculation->state_data();
  auto const &state = *state_data.state;
  auto const &composition_calculator = get_composition_calculator(system);
  auto const &orbits =
      get_local_basis_set_cluster_info(system, m_data->local_basis_set_name)
          ->orbits;
  auto prim_nlist = system.prim_neighbor_list;
  auto supercell_nlist = get_supercell_neighbor_list(system, state);
  auto const &supercell_index_converter =
      get_index_conversions(system, state).index_converter();
  clexulator::ConfigDoFValues const *dof_values = &get_dof_values(state);

  m_local_orbit_composition_calculator =
      std::make_shared<LocalOrbitCompositionCalculator>(
          orbits, m_data->orbits_to_calculate, m_data->combine_orbits,
          prim_nlist, supercell_nlist, supercell_index_converter,
          composition_calculator, dof_values);
}

/// \brief Make local orbit composition collecting functions
/// ("local_orbit_composition.<key>")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_local_orbit_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  auto const &system = get_system(calculation);
  std::vector<monte::DiscreteVectorIntHistogramFunction> f_list;

  for (auto const &pair : system.local_orbit_composition_calculator_data) {
    std::string key = pair.first;
    std::shared_ptr<LocalOrbitCompositionCalculatorData const> data =
        pair.second;
    std::string event_type_name = data->event_type_name;
    if (event_type_name.empty()) {
      continue;
    }

    SelectedEventInfo info = make_selected_event_info(calculation);
    info.make_indices_by_equivalent_index_per_event_type(event_type_name);
    auto prim_event_index_to_index = info.prim_event_index_to_index;
    auto prim_event_index_to_has_value = info.prim_event_index_to_has_value;
    auto selected_event = calculation->selected_event();

    LocalOrbitCompositionCollector collector(calculation, data);

    jsonParser tjson;
    std::stringstream desc;
    desc << "Selected event local orbit composition calculator " << key
         << " for event=" << data->event_type_name
         << ", local_basis_set=" << data->local_basis_set_name
         << ", orbits=" << to_json(data->orbits_to_calculate, tjson)
         << ", combine_orbits=" << std::boolalpha << data->combine_orbits;

    // Create the function
    monte::DiscreteVectorIntHistogramFunction f(
        /* name */
        "local_orbit_composition." + key,
        /* description */
        desc.str(),
        /* shape */
        collector.shape(),
        /* component_names */
        collector.component_names(),
        /* requires_event_state */
        false,
        /* function */
        [selected_event, collector]() mutable /*allow collector to change*/
        -> Eigen::VectorXl {
          return collector
              .value(get_unitcell_index(selected_event),
                     get_equivalent_index(selected_event))
              .reshaped()
              .cast<long>();
        },
        /* has_value_function */
        [prim_event_index_to_has_value, selected_event]() -> bool {
          return (*prim_event_index_to_has_value)[get_prim_event_index(
              selected_event)];
        },
        /* max_size */
        data->max_size);
    f_list.push_back(f);
  }

  return f_list;
}

/// \brief Make dE_activated collecting function
/// ("dE_activated.by_type")
monte::PartitionedHistogramFunction<double> make_dE_activated_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  SelectedEventInfo info = make_selected_event_info(calculation);
  info.make_indices_by_type();
  auto prim_event_index_to_index = info.prim_event_index_to_index;
  auto selected_event = calculation->selected_event();

  // Create the function
  monte::PartitionedHistogramFunction<double> f(
      /* name */
      "dE_activated.by_type",
      /* description */
      "Selected event activated state energy, relative to the initial state, "
      "partitioned by event type",
      /* requires_event_state */
      true,
      /* function */
      [selected_event]() {
        return get_event_state(selected_event).dE_activated;
      },
      /* partition_names */
      info.partition_names,
      /* get_partition */
      [prim_event_index_to_index, selected_event]() {
        return (
            *prim_event_index_to_index)[get_prim_event_index(selected_event)];
      },
      /* is_log */
      false,
      /* initial_begin */
      0.0,
      /* bin_width */
      0.01,
      /* max_size */
      1000);

  return f;
}

/// \brief Make dE_activated collecting function
/// ("dE_activated.by_equivalent_index")
monte::PartitionedHistogramFunction<double>
make_dE_activated_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  SelectedEventInfo info = make_selected_event_info(calculation);
  info.make_indices_by_equivalent_index();
  auto prim_event_index_to_index = info.prim_event_index_to_index;
  auto selected_event = calculation->selected_event();

  // Create the function
  monte::PartitionedHistogramFunction<double> f(
      /* name */
      "dE_activated.by_equivalent_index",
      /* description */
      "Selected event activated state energy, relative to the initial state, "
      "partitioned by event type and equivalent index.",
      /* requires_event_state */
      true,
      /* function */
      [selected_event]() {
        return get_event_state(selected_event).dE_activated;
      },
      /* partition_names */
      info.partition_names,
      /* get_partition */
      [prim_event_index_to_index, selected_event]() {
        return (
            *prim_event_index_to_index)[get_prim_event_index(selected_event)];
      },
      /* is_log */
      false,
      /* initial_begin */
      0.0,
      /* bin_width */
      0.01,
      /* max_size */
      1000);

  return f;
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM
