#ifndef CASM_clexmonte_monte_calculator_KineticCalculator
#define CASM_clexmonte_monte_calculator_KineticCalculator

#include "casm/clexmonte/kinetic/kinetic_events.hh"
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/monte/MethodLog.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic_2 {

class KineticPotential : public BaseMontePotential {
 public:
  KineticPotential(std::shared_ptr<StateData> _state_data);

  // --- Data used in the potential calculation: ---

  state_type const &state;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() override;

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() override;

  /// \brief Calculate change in (per_supercell) potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) override;
};

enum class kinetic_event_data_type {
  high_memory,
  default_memory,
  low_memory,
};

class KineticCalculator : public BaseMonteCalculator {
 public:
  using BaseMonteCalculator::engine_type;

  KineticCalculator();

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override;

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override;

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override;

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override;

  /// \brief Construct functions that may be used to collect selected event data
  std::optional<monte::SelectedEventDataFunctions>
  standard_selected_event_data_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override;

  /// \brief Construct default SamplingFixtureParams
  sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results, bool write_trajectory, bool write_observations,
      bool write_status, std::optional<std::string> output_dir,
      std::optional<std::string> log_file,
      double log_frequency_in_s) const override;

  /// \brief Validate the state's configuration
  Validator validate_configuration(state_type &state) const override;

  /// \brief Validate state's conditions
  Validator validate_conditions(state_type &state) const override;

  /// \brief Validate state
  Validator validate_state(state_type &state) const override;

  /// \brief Validate and set the current state, construct state_data, construct
  ///     potential
  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) override;

  /// \brief Set event data (includes calculating all rates), using current
  /// state data
  void set_event_data(std::shared_ptr<engine_type> engine) override;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override;

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           run_manager_type<engine_type> &run_manager) override;

  // --- KineticCalculator specific functions ---

  /// \brief Print a warning to std::cerr if events with no barrier were
  ///     encountered
  void check_n_not_normal(
      std::map<std::string, Index> const &n_not_normal) const;

  // --- Parameters ---

  // Verbosity level (if applicable)
  int verbosity_level = 10;

  // Used by state modifying functions
  double mol_composition_tol = CASM::TOL;

  // Type of event data structure
  kinetic_event_data_type event_data_type =
      kinetic_event_data_type::default_memory;

  // Type of event selector
  kinetic_event_selector_type event_selector_type =
      kinetic_event_selector_type::vector_sum_tree;

  // Type of impact table:
  // - Only takes effect if event_data_type is `default_memory` or `low_memory`
  // - If true: somewhat higher memory use; somewhat faster impact list
  // - If false: somewhat lower memory use; somewhat slower impact list
  bool use_neighborlist_impact_table = true;

  // If true, events without barriers are allowed with warning messages;
  // If false (default), an exception is thrown at the `select_event` step of
  // a run if an event without a barrier is encountered
  bool allow_events_with_no_barrier = false;

  /// If true (default) check if potentially impacted events are allowed
  /// and only assign them to the event list if they are (adds an
  /// additional check, but may reduce the size of the event list).
  /// Otherwise, assign all potentially impacted events to the event list
  /// (whether they are allowed will still be checked during the rate
  /// calculation).
  bool assign_allowed_events_only = true;

  /// \brief Reset the derived Monte Carlo calculator
  void _reset() override;

  /// \brief Clone the KineticCalculator
  KineticCalculator *_clone() const override;
};

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator();
}

#endif