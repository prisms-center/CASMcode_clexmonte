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
  low_memory, /* currently not used */
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
  std::optional<monte::SelectedEventFunctions>
  standard_selected_event_functions(
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
  void set_event_data() override;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override;

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           run_manager_type<engine_type> &run_manager) override;

  // --- KineticCalculator specific functions ---

  /// \brief Print a warning to std::cerr if abnormal events were
  ///     encountered
  void check_n_encountered_abnormal(
      std::map<std::string, Index> const &n_encountered_abnormal) const;

  /// \brief Print a warning to std::cerr if abnormal events were
  ///     selected
  void check_n_selected_abnormal(
      std::map<std::string, Index> const &n_selected_abnormal) const;

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

  // Event filters
  std::optional<std::vector<EventFilterGroup>> event_filters;

  // Print event data summary (default=false)
  bool print_event_data_summary = false;

  // Various options for event list implementation and event handling
  EventDataOptions event_data_options;

  /// \brief Reset the derived Monte Carlo calculator
  void _reset() override;

  /// \brief Clone the KineticCalculator
  KineticCalculator *_clone() const override;

  template <bool DebugMode>
  void make_complete_event_data_impl();

  template <bool DebugMode>
  void make_allowed_event_data_impl();
};

}  // namespace kinetic_2
}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator();
}

#endif