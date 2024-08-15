#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/events/lotto/rejection_free.hpp"
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/analysis_functions.hh"
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"
#include "casm/clexmonte/run/functions.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/monte/methods/nfold.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

namespace CASM {

namespace clexmonte {

occ_events::OccPosition _make_atom_position(
    occ_events::OccSystem const &event_system,
    xtal::UnitCellCoord const &integral_site_coordinate, Index occupant_index) {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= event_system.prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "integral_site_coordinate");
  }
  std::vector<xtal::Molecule> const &occupant_dof =
      event_system.prim->basis()[b].occupant_dof();
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid occupant_index");
  }

  return occ_events::OccPosition{false, true, integral_site_coordinate,
                                 occupant_index, 0};
}

/// \brief Convert monte::OccCanditate to occ_events::OccPosition
occ_events::OccPosition _make_atom_position(
    occ_events::OccSystem const &event_system,
    monte::Conversions const &convert, monte::OccCandidate const &cand) {
  Index b = *convert.asym_to_b(cand.asym).begin();
  return _make_atom_position(event_system, xtal::UnitCellCoord(b, 0, 0, 0),
                             convert.occ_index(cand.asym, cand.species_index));
}

/// \brief Make atom-in-resevoir position for same chemical type as given
/// OccPosition
occ_events::OccPosition _make_atom_in_resevoir_position(
    occ_events::OccSystem const &event_system,
    occ_events::OccPosition const &pos) {
  Index chemical_index = event_system.get_chemical_index(pos);
  return occ_events::OccPosition{true, true, xtal::UnitCellCoord{0, 0, 0, 0},
                                 chemical_index, 0};
}

/// \brief Convert monte::OccSwap to occ_events::OccEvent
occ_events::OccEvent _make_semigrand_canonical_swap_event(
    monte::OccSwap const &swap, occ_events::OccSystem const &event_system,
    monte::Conversions const &convert) {
  occ_events::OccPosition pos_a =
      _make_atom_position(event_system, convert, swap.cand_a);
  occ_events::OccPosition resevoir_a =
      _make_atom_in_resevoir_position(event_system, pos_a);

  occ_events::OccPosition pos_b =
      _make_atom_position(event_system, convert, swap.cand_b);
  occ_events::OccPosition resevoir_b =
      _make_atom_in_resevoir_position(event_system, pos_b);

  occ_events::OccTrajectory traj_a({pos_a, resevoir_a});
  occ_events::OccTrajectory traj_b({resevoir_b, pos_b});

  return occ_events::OccEvent({traj_a, traj_b});
}
/// \brief Make OccEventTypeData for semi-grand canonical events
std::map<std::string, OccEventTypeData> _make_event_type_data(
    std::shared_ptr<system_type> system, state_type const &state,
    std::vector<monte::OccSwap> const &semigrand_canonical_swaps) {
  auto event_system = get_event_system(*system);
  monte::Conversions const &convert = get_index_conversions(*system, state);

  std::map<std::string, OccEventTypeData> event_type_data;
  auto const &occevent_symgroup_rep = get_occevent_symgroup_rep(*system);
  for (monte::OccSwap const &swap : semigrand_canonical_swaps) {
    // do not repeat forward and reverse -
    //   reverse will be constructed by make_prim_event_list(
    if (swap.cand_a.species_index < swap.cand_b.species_index) {
      continue;
    }
    occ_events::OccEvent event =
        _make_semigrand_canonical_swap_event(swap, *event_system, convert);
    std::set<occ_events::OccEvent> orbit =
        occ_events::make_prim_periodic_orbit(event, occevent_symgroup_rep);

    std::string event_type_name =
        "swap-" + std::to_string(swap.cand_a.asym) + "-" +
        std::to_string(swap.cand_a.species_index) + "-" +
        std::to_string(swap.cand_b.species_index);
    event_type_data[event_type_name].events =
        std::vector<occ_events::OccEvent>(orbit.begin(), orbit.end());
  }
  return event_type_data;
}

// The NfoldPotential is the same as the SemiGrandCanonicalPotential
class NfoldPotential : public BaseMontePotential {
 public:
  NfoldPotential(std::shared_ptr<StateData> _state_data)
      : BaseMontePotential(_state_data),
        state(*state_data->state),
        n_unitcells(state_data->n_unitcells),
        occupation(get_occupation(state)),
        convert(*state_data->convert),
        composition_calculator(
            get_composition_calculator(*this->state_data->system)),
        composition_converter(
            get_composition_converter(*this->state_data->system)),
        param_chem_pot(state.conditions.vector_values.at("param_chem_pot")),
        formation_energy_clex(
            get_clex(*state_data->system, state, "formation_energy")) {
    if (param_chem_pot.size() !=
        composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalPotential: param_chem_pot size error");
    }

    exchange_chem_pot =
        make_exchange_chemical_potential(param_chem_pot, composition_converter);
  }

  // --- Data used in the potential calculation: ---

  state_type const &state;
  Index n_unitcells;
  Eigen::VectorXi const &occupation;
  monte::Conversions const &convert;
  composition::CompositionCalculator const &composition_calculator;
  composition::CompositionConverter const &composition_converter;
  Eigen::VectorXd param_chem_pot;
  std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex;
  Eigen::MatrixXd exchange_chem_pot;

  /// \brief Calculate (per_supercell) potential value
  double per_supercell() override {
    Eigen::VectorXd mol_composition =
        composition_calculator.mean_num_each_component(occupation);
    Eigen::VectorXd param_composition =
        composition_converter.param_composition(mol_composition);

    return formation_energy_clex->per_supercell() -
           n_unitcells * param_chem_pot.dot(param_composition);
  }

  /// \brief Calculate (per_unitcell) potential value
  double per_unitcell() override { return this->per_supercell() / n_unitcells; }

  /// \brief Calculate change in (per_supercell) semi-grand potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ) override {
    double delta_formation_energy =
        formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
    double delta_potential_energy = delta_formation_energy;
    for (Index i = 0; i < linear_site_index.size(); ++i) {
      Index l = linear_site_index[i];
      Index asym = convert.l_to_asym(l);
      Index curr_species = convert.species_index(asym, occupation(l));
      Index new_species = convert.species_index(asym, new_occ[i]);
      delta_potential_energy -= exchange_chem_pot(new_species, curr_species);
    }

    return delta_potential_energy;
  }
};

/// NfoldEventData data --> These are Nfold specific and ported over from
/// include/casm/clexmonte/nfold/nfold_events.hh
/// \brief Data calculated for a single event in a single state
struct EventState {
  bool is_allowed;  ///< Is allowed given current configuration
  double dE_final;  ///< Final state energy, relative to initial state
  double rate;      ///< Occurance rate
};

/// \brief CompleteEventCalculator is an event calculator with the required
///     interface for the classes `lotto::RejectionFree` and `lotto::Rejection`.
///
/// Notes:
/// - Expected to be constructed as shared_ptr
/// - Mostly holds references to external data structures
/// - Stores one `EventState` which is used to perform the calculations
struct CompleteEventCalculator {
  CompleteEventCalculator(std::shared_ptr<BaseMontePotential> _potential,
                          std::vector<PrimEventData> const &_prim_event_list,
                          std::map<EventID, EventData> const &_event_list)
      : prim_event_list(_prim_event_list),
        event_list(_event_list),
        potential(_potential) {}
  /// \brief Prim event list
  std::vector<PrimEventData> const &prim_event_list;

  /// \brief Complete event list
  std::map<EventID, EventData> const &event_list;

  /// \brief Holds last calculated event state
  EventState event_state;

  /// \brief Potential
  std::shared_ptr<BaseMontePotential> potential;

  /// \brief Get CASM::monte::OccEvent corresponding to given event ID
  double calculate_rate(EventID const &id) {
    EventData const &event_data = event_list.at(id);
    PrimEventData const &prim_event_data =
        prim_event_list.at(id.prim_event_index);

    /// ---
    /// CHECK: Does getting dof_values make sense?
    clexulator::ConfigDoFValues const *dof_values =
        &potential->state_data->state->configuration.dof_values;

    int i = 0;
    for (Index l : event_data.event.linear_site_index) {
      if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
        event_state.is_allowed = false;
        event_state.rate = 0.0;
        return event_state.rate;
      }
      ++i;
    }
    event_state.is_allowed = true;

    // calculate change in energy to final state
    event_state.dE_final = potential->occ_delta_per_supercell(
        event_data.event.linear_site_index, prim_event_data.occ_final);

    // calculate rate
    // CHECK: does getting beta make sense?
    if (event_state.dE_final <= 0.0) {
      event_state.rate = 1.0;
    } else {
      event_state.rate = exp(
          -potential->state_data->state->conditions.scalar_values.at("beta") *
          event_state.dE_final);
    }

    /// ---
    return event_state.rate;
  }
};

struct NfoldEventData {
  NfoldEventData(std::shared_ptr<system_type> system, state_type const &state,
                 monte::OccLocation const &occ_location,
                 std::vector<monte::OccSwap> const &semigrand_canonical_swaps,
                 std::shared_ptr<BaseMontePotential> potential) {
    // Make OccEvents from SemiGrandCanonical swaps
    // key: event_type_name, value: symmetrically equivalent events
    system->event_type_data =
        _make_event_type_data(system, state, semigrand_canonical_swaps);

    prim_event_list = clexmonte::make_prim_event_list(*system);

    prim_impact_info_list = clexmonte::make_prim_impact_info_list(
        *system, prim_event_list, {"formation_energy"});

    // TODO: rejection-clexmonte option does not require impact table
    event_list = clexmonte::make_complete_event_list(
        prim_event_list, prim_impact_info_list, occ_location);

    // Construct CompleteEventCalculator
    event_calculator = std::make_shared<CompleteEventCalculator>(
        potential, prim_event_list, event_list.events);
  }

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> prim_event_list;

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  /// All supercell events, and which events must be updated
  /// when one occurs
  clexmonte::CompleteEventList event_list;

  /// Calculator for KMC event selection
  std::shared_ptr<CompleteEventCalculator> event_calculator;
};

// --------------------------------------------------------------------------
// NfoldCalculator and implementation----------------------------------------
class NfoldCalculator : public BaseMonteCalculator {
 public:
  using BaseMonteCalculator::engine_type;

  NfoldCalculator()
      : BaseMonteCalculator("NfoldCalculator",     // calculator_name
                            {},                    // required_basis_set,
                            {},                    // required_local_basis_set,
                            {"formation_energy"},  // required_clex,
                            {},                    // required_multiclex,
                            {},                    // required_local_clex,
                            {},                    // required_local_multiclex,
                            {},                    // required_dof_spaces,
                            {},                    // required_params,
                            {},                    // optional_params,
                            true,                  // time_sampling_allowed,
                            false,                 // update_species,
                            false                  // is_multistate_method,
        ) {}

  /// Hold event_data. This needs to be in BaseMonteCalculator.hh to be exposed
  /// to the python side
  std::shared_ptr<NfoldEventData> event_data;

  /// Hold NfoldData
  monte::NfoldData<config_type, statistics_type, engine_type> nfold_data;

  // The sampling_functions, json_sampling_functions, analysis_functions,
  // modifying_functions are all copied over from
  // SemiGrandCanonicalMonteCalculator
  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, state_sampling_function_type>
  standard_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<state_sampling_function_type> functions =
        monte_calculator::common_sampling_functions(
            calculation, "potential_energy",
            "Potential energy of the state (normalized per primitive cell)");

    // Specific to semi-grand canonical
    functions.push_back(monte_calculator::make_param_chem_pot_f(calculation));

    std::map<std::string, state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to sample various quantities
  ///     of the Monte Carlo calculation as it runs
  std::map<std::string, json_state_sampling_function_type>
  standard_json_sampling_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<json_state_sampling_function_type> functions =
        monte_calculator::common_json_sampling_functions(calculation);

    std::map<std::string, json_state_sampling_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to analyze Monte Carlo
  ///     calculation results
  std::map<std::string, results_analysis_function_type>
  standard_analysis_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    std::vector<results_analysis_function_type> functions = {
        monte_calculator::make_heat_capacity_f(calculation),
        monte_calculator::make_mol_susc_f(calculation),
        monte_calculator::make_param_susc_f(calculation),
        monte_calculator::make_mol_thermochem_susc_f(calculation),
        monte_calculator::make_param_thermochem_susc_f(calculation)};

    std::map<std::string, results_analysis_function_type> function_map;
    for (auto const &f : functions) {
      function_map.emplace(f.name, f);
    }
    return function_map;
  }

  /// \brief Construct functions that may be used to modify states
  StateModifyingFunctionMap standard_modifying_functions(
      std::shared_ptr<MonteCalculator> const &calculation) const override {
    return StateModifyingFunctionMap();
  }

  /// \brief Construct default SamplingFixtureParams
  sampling_fixture_params_type make_default_sampling_fixture_params(
      std::shared_ptr<MonteCalculator> const &calculation, std::string label,
      bool write_results, bool write_trajectory, bool write_observations,
      bool write_status, std::optional<std::string> output_dir,
      std::optional<std::string> log_file,
      double log_frequency_in_s) const override {
    monte::SamplingParams sampling_params;
    {
      auto &s = sampling_params;
      s.sampler_names = {"clex.formation_energy", "potential_energy",
                         "mol_composition", "param_composition"};
      std::string prefix;
      prefix = "order_parameter_";
      for (auto const &pair : calculation->system()->dof_spaces) {
        s.sampler_names.push_back(prefix + pair.first);
      }
      prefix = "subspace_order_parameter_";
      for (auto const &pair : calculation->system()->dof_subspaces) {
        s.sampler_names.push_back(prefix + pair.first);
      }
      if (write_trajectory) {
        s.do_sample_trajectory = true;
      }
    }

    monte::CompletionCheckParams<statistics_type> completion_check_params;
    {
      auto &c = completion_check_params;
      c.equilibration_check_f = monte::default_equilibration_check;
      c.calc_statistics_f =
          monte::default_statistics_calculator<statistics_type>();

      converge(calculation->sampling_functions, completion_check_params)
          .set_abs_precision("potential_energy", 0.001)
          .set_abs_precision("param_composition", 0.001);
    }

    std::vector<std::string> analysis_names = {
        "heat_capacity", "mol_susc", "param_susc", "mol_thermochem_susc",
        "param_thermochem_susc"};

    return clexmonte::make_sampling_fixture_params(
        label, calculation->sampling_functions,
        calculation->json_sampling_functions, calculation->analysis_functions,
        sampling_params, completion_check_params, analysis_names, write_results,
        write_trajectory, write_observations, write_status, output_dir,
        log_file, log_frequency_in_s);
  }

  /// \brief Validate the state's configuration (all are valid)
  Validator validate_configuration(state_type &state) const override {
    return Validator{};
  }

  /// \brief Validate state's conditions
  ///
  /// Notes:
  /// - requires scalar temperature
  /// - requires vector param_chem_pot
  /// - warnings if other conditions are present
  Validator validate_conditions(state_type &state) const override {
    // validate state.conditions
    monte::ValueMap const &conditions = state.conditions;
    Validator v;
    v.insert(validate_keys(conditions.scalar_values,
                           {"temperature"} /*required*/, {} /*optional*/,
                           "scalar", "condition", false /*throw_if_invalid*/));
    v.insert(validate_keys(conditions.vector_values,
                           {"param_chem_pot"} /*required*/, {} /*optional*/,
                           "vector", "condition", false /*throw_if_invalid*/));

    return v;
  }

  /// \brief Validate state
  Validator validate_state(state_type &state) const override {
    Validator v;
    v.insert(this->validate_configuration(state));
    v.insert(this->validate_conditions(state));
    return v;
  }

  void set_state_and_potential(state_type &state,
                               monte::OccLocation *occ_location) override {
    // Validate system
    if (this->system == nullptr) {
      throw std::runtime_error(
          "Error in NfoldCalculator::run: system==nullptr");
    }

    // Validate state
    Validator v = this->validate_state(state);
    print(CASM::log(), v);
    if (!v.valid()) {
      throw std::runtime_error(
          "Error in NfoldCalculator::run: Invalid initial state");
    }

    // Make state data
    this->state_data =
        std::make_shared<StateData>(this->system, &state, occ_location);

    // Make potential calculator
    this->potential = std::make_shared<NfoldPotential>(this->state_data);
  }

  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<engine_type> &run_manager) override {
    this->set_state_and_potential(state, &occ_location);

    // Get semi-grand canonical swaps
    std::vector<monte::OccSwap> const &semigrand_canonical_swaps =
        get_semigrand_canonical_swaps(*this->system);

    // CHECK: In the older implementation, it does not make the event list
    // if supercell of the state is the same as the supercell of the calculator
    // How does that work here?
    this->event_data = std::make_shared<NfoldEventData>(
        this->system, state, occ_location, semigrand_canonical_swaps,
        this->potential);

    // nfold data
    monte::Conversions const &convert =
        get_index_conversions(*this->system, state);
    Index n_allowed_per_unitcell =
        get_n_allowed_per_unitcell(convert, semigrand_canonical_swaps);
    this->nfold_data.n_events_possible =
        static_cast<double>(this->state_data->n_unitcells) *
        n_allowed_per_unitcell;

    // Make selector
    lotto::RejectionFreeEventSelector event_selector(
        this->event_data->event_calculator,
        clexmonte::make_complete_event_id_list(
            this->state_data->n_unitcells, this->event_data->prim_event_list),
        this->event_data->event_list.impact_table,
        std::make_shared<lotto::RandomGenerator>(run_manager.engine));

    auto get_event_f = [&](EventID const &selected_event_id) {
      return this->event_data->event_list.events.at(selected_event_id).event;
    };

    // Run nfold-way
    monte::nfold<EventID>(state, occ_location, this->nfold_data, event_selector,
                          get_event_f, run_manager);
  }

  /// \brief Perform a single run, evolving one or more states
  void run(int current_state, std::vector<state_type> &states,
           std::vector<monte::OccLocation> &occ_locations,
           run_manager_type<engine_type> &run_manager) override {
    throw std::runtime_error(
        "Error: SemiGrandCanonicalCalculator does not allow multi-state runs");
  }

  // --- Parameters ---
  int verbosity_level = 10;

  /// \brief Reset the derived Monte Carlo calculator
  ///
  /// Parameters:
  ///
  ///   verbosity: str or int, default=10
  ///       If integer, the allowed range is `[0,100]`. If string, then:
  ///       - "none" is equivalent to integer value 0
  ///       - "quiet" is equivalent to integer value 5
  ///       - "standard" is equivalent to integer value 10
  ///       - "verbose" is equivalent to integer value 20
  ///       - "debug" is equivalent to integer value 100
  void _reset() override {
    ParentInputParser parser{params};

    // "verbosity": str or int, default=10
    this->verbosity_level = parse_verbosity(parser);
    CASM::log().set_verbosity(this->verbosity_level);

    // TODO: enumeration

    std::stringstream ss;
    ss << "Error in SemiGrandCanonicalCalculator: error reading calculation "
          "parameters.";
    std::runtime_error error_if_invalid{ss.str()};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    return;
  }

  /// \brief Clone the SemiGrandCanonicalCalculator
  NfoldCalculator *_clone() const override {
    return new NfoldCalculator(*this);
  }
};

}  // namespace clexmonte
}  // namespace CASM

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// NfoldCalculator
CASM::clexmonte::BaseMonteCalculator *make_NfoldCalculator() {
  return new CASM::clexmonte::NfoldCalculator();
}
}
