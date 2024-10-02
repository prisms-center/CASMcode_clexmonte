/// An implementation of a kinetic Monte Carlo main loop
/// that makes use of the RunManager provided by
/// casm/monte/run_management to implement sampling
/// fixtures and results data structures and input/output
/// methods and a data structure that allows sampling
/// atomic displacements for kinetic coefficient
/// calculations.

#ifndef CASM_clexmonte_methods_kinetic_monte_carlo
#define CASM_clexmonte_methods_kinetic_monte_carlo

// logging
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/MethodLog.hh"
#include "casm/monte/checks/io/json/CompletionCheck_json_io.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/methods/kinetic_monte_carlo.hh"
#include "casm/monte/run_management/State.hh"
#include "casm/monte/sampling/SelectedEventData.hh"

namespace CASM {
namespace clexmonte {

/// \brief Construct a list of atom names corresponding to OccLocation atoms
std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system);

template <typename EventIDType, typename ConfigType,
          typename SetSelectedEventFunction, typename StatisticsType,
          typename EngineType>
void kinetic_monte_carlo_v2(
    state_type &state, monte::OccLocation &occ_location,
    monte::KMCData<ConfigType, StatisticsType, EngineType> &kmc_data,
    SelectedEvent &selected_event,
    SetSelectedEventFunction set_selected_event_f,
    std::optional<monte::SelectedEventDataCollector> &collector,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system);

// --- Implementation ---

/// \brief Construct a list of atom names corresponding to OccLocation atoms
///
/// Notes:
/// - If atoms are conserved, then the order of this list will remain unchanged
///   during the course of a calculation
/// - Values are set to -1 if atom is no longer in supercell
inline std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system) {
  // sanity check:
  monte::Conversions const &convert = occ_location.convert();
  if (convert.species_size() != occ_system.orientation_name_list.size()) {
    throw std::runtime_error(
        "Error in CASM::clexmonte::kinetic::make_snapshot_for_conserved_atoms: "
        "mismatch between monte::Conversions and occ_events::OccSystem.");
  }

  // collect atom name indices
  std::vector<Index> atom_name_index_list(occ_location.atom_size(), -1);
  for (Index i = 0; i < occ_location.mol_size(); ++i) {
    monte::Mol const &mol = occ_location.mol(i);
    Index b = convert.l_to_b(mol.l);
    Index occupant_index =
        occ_system.orientation_to_occupant_index[b][mol.species_index];
    Index atom_position_index = 0;
    for (Index atom_id : mol.component) {
      Index atom_name_index =
          occ_system.atom_position_to_name_index[b][occupant_index]
                                                [atom_position_index];
      atom_name_index_list.at(atom_id) = atom_name_index;
      ++atom_position_index;
    }
  }
  return atom_name_index_list;
}

/// \brief Run a kinetic Monte Carlo calculation
///
/// \param state The state. Consists of both the initial
///     configuration and conditions. Conditions must include `temperature`
///     and any others required by `potential`.
/// \param occ_location An occupant location tracker, which enables efficient
///     event proposal. It must already be initialized with the input state.
/// \param kmc_data Stores data to be made available to the sampling functions
///     along with the current state.
/// \param selected_event The selected event data object which the selected
///     event data functions are expecting to be set with the selected event
///     data.
/// \param set_selected_event_f A function that can set `selected_event` with
///     signature `set_selected_event_f(SelectedEvent &selected_event, bool
///     requires_event_state)`.
/// \param collector Collects selected event data
/// \param run_manager Contains sampling fixtures and after completion holds
///     final results
/// \param event_system Defines the system for OccPosition / OccTrajectory /
///     OccEvent. Used in particular to determine
///     `kmc_data.atom_name_index_list`.
///
/// \returns A Results<ConfigType> instance with run results.
///
/// Required interface for `State<ConfigType>`:
/// - `Eigen::VectorXi &get_occupation(State<ConfigType> const &state)`
/// - `Eigen::Matrix3l const &get_transformation_matrix_to_super(
///        State<ConfigType> const &state)`
///
/// State properties that are set:
/// - None
///
template <typename EventIDType, typename ConfigType,
          typename SetSelectedEventFunction, typename StatisticsType,
          typename EngineType>
void kinetic_monte_carlo_v2(
    state_type &state, monte::OccLocation &occ_location,
    monte::KMCData<ConfigType, StatisticsType, EngineType> &kmc_data,
    SelectedEvent &selected_event,
    SetSelectedEventFunction set_selected_event_f,
    std::optional<monte::SelectedEventDataCollector> &collector,
    monte::RunManager<ConfigType, StatisticsType, EngineType> &run_manager,
    std::shared_ptr<occ_events::OccSystem> event_system) {
  // Validate existence of sampling fixtures
  if (run_manager.sampling_fixtures.empty()) {
    throw std::runtime_error(
        "Error in clexmonte::kinetic_monte_carlo_v2: no sampling fixtures");
  }
  monte::MonteCounter const &counter =
      run_manager.sampling_fixtures.front()->counter();

  // Initialize atom positions & time
  kmc_data.sampling_fixture_label.clear();
  kmc_data.sampling_fixture = nullptr;
  kmc_data.total_rate = 0.0;
  kmc_data.time = 0.0;
  kmc_data.unique_atom_id = occ_location.unique_atom_id();
  kmc_data.atom_positions_cart = occ_location.atom_positions_cart();
  kmc_data.prev_time.clear();
  kmc_data.prev_atom_positions_cart.clear();
  kmc_data.prev_unique_atom_id.clear();
  for (auto &fixture_ptr : run_manager.sampling_fixtures) {
    kmc_data.prev_time.emplace(fixture_ptr->label(), kmc_data.time);
    kmc_data.prev_atom_positions_cart.emplace(fixture_ptr->label(),
                                              kmc_data.atom_positions_cart);
    kmc_data.prev_unique_atom_id.emplace(fixture_ptr->label(),
                                         kmc_data.unique_atom_id);
  }

  // Pre- and post- sampling actions:

  // notes: it is important this uses
  // - the total_rate obtained for the state before applying the selected event
  auto pre_sample_action =
      [&](monte::SamplingFixture<ConfigType, StatisticsType, EngineType>
              &fixture,
          monte::State<ConfigType> const &state) {
        // set data that can be used in sampling functions
        kmc_data.sampling_fixture_label = fixture.label();
        kmc_data.sampling_fixture = &fixture;
        kmc_data.unique_atom_id = occ_location.unique_atom_id();
        kmc_data.atom_name_index_list =
            make_atom_name_index_list(occ_location, *event_system);
        kmc_data.atom_positions_cart = occ_location.atom_positions_cart();
        kmc_data.total_rate = selected_event.total_rate;
        if (fixture.params().sampling_params.sample_mode ==
            monte::SAMPLE_MODE::BY_TIME) {
          kmc_data.time = fixture.next_sample_time();
        }
      };

  auto post_sample_action =
      [&](monte::SamplingFixture<ConfigType, StatisticsType, EngineType>
              &fixture,
          monte::State<ConfigType> const &state) {
        // set data that can be used in sampling functions
        kmc_data.prev_time[fixture.label()] = kmc_data.time;
        kmc_data.prev_atom_positions_cart[fixture.label()] =
            kmc_data.atom_positions_cart;
        kmc_data.prev_unique_atom_id[fixture.label()] = kmc_data.unique_atom_id;
      };

  // Main loop
  double event_time;
  bool collect_selected_event_data = collector.has_value();
  bool requires_event_state =
      collector.has_value() && collector->requires_event_state;
  selected_event.reset();
  run_manager.initialize(occ_location.mol_size());
  run_manager.update_next_sampling_fixture();
  while (!run_manager.is_complete()) {
    run_manager.write_status_if_due();

    // Select an event. This function:
    // - Updates rates of events impacted by the previous selected event (if
    //   there was a previous event)
    // - Updates the total rate
    // - Chooses an event and time increment (does not apply event)
    // - Sets a list of impacted events by the chosen event that will be updated
    //   on the next iteration
    // - If `requires_event_state` is true, then the event state is calculated
    //   for the selected event
    set_selected_event_f(selected_event, requires_event_state);
    event_time = kmc_data.time + selected_event.time_increment;

    // Sample data, if a sample is due by count
    // - This location correctly handles sampling at count=0 and count!=0
    // - If the sample count is n, then the state after the n-th step/pass is
    //   sampled.
    run_manager.sample_data_by_count_if_due(state, pre_sample_action,
                                            post_sample_action);

    // Sample data, if a sample is due by time
    // - This location correctly handles sampling at sample times >= 0
    // - Samples are taken using the state before the event for each
    //   sample time <= the event time
    // - If the sample time is exactly equal to the event time (should be
    //   vanishingly rare), then the state before the event occurs is sampled.
    run_manager.sample_data_by_time_if_due(event_time, state, pre_sample_action,
                                           post_sample_action);

    // Set time -- for all fixtures
    run_manager.set_time(event_time);

    // Increment count -- for all fixtures
    run_manager.increment_step();

    // Collect selected event data
    if (collect_selected_event_data) {
      collector->collect(counter);
    }

    // Apply event
    run_manager.increment_n_accept();
    occ_location.apply(selected_event.event_data->event, get_occupation(state));
    kmc_data.time = event_time;
  }

  run_manager.finalize(state);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
