#ifndef CASM_clexmonte_nfold_impl
#define CASM_clexmonte_nfold_impl

#include "casm/clexmonte/nfold/nfold.hh"
#include "casm/clexmonte/semigrand_canonical/calculator_impl.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/nfold.hh"

namespace CASM {
namespace clexmonte {
namespace nfold {

template <typename EngineType>
Nfold<EngineType>::Nfold(std::shared_ptr<system_type> _system)
    : semigrand_canonical::SemiGrandCanonical<EngineType>(_system) {}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void Nfold<EngineType>::run(state_type &state, monte::OccLocation &occ_location,
                            run_manager_type<EngineType> &run_manager) {
  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("param_chem_pot")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `param_chem_pot` conditions not set.");
  }

  // Store state info / pointers
  this->state = &state;
  this->occ_location = &occ_location;
  this->conditions =
      std::make_shared<semigrand_canonical::SemiGrandCanonicalConditions>(
          get_composition_converter(*this->system));
  this->conditions->set_all(state.conditions, false);
  Index n_unitcells = this->transformation_matrix_to_super.determinant();

  // Make potential calculator
  this->potential =
      std::make_shared<semigrand_canonical::SemiGrandCanonicalPotential>(
          this->system);
  this->potential->set(this->state, this->conditions);
  this->formation_energy = this->potential->formation_energy();

  // Get swaps
  std::vector<monte::OccSwap> const &semigrand_canonical_swaps =
      get_semigrand_canonical_swaps(*this->system);

  // if same supercell
  // -> just re-set potential & avoid re-constructing event list
  if (this->transformation_matrix_to_super ==
          get_transformation_matrix_to_super(state) &&
      this->conditions != nullptr) {
    this->event_data->event_calculator->potential = this->potential;
  } else {
    this->transformation_matrix_to_super =
        get_transformation_matrix_to_super(state);
    n_unitcells = this->transformation_matrix_to_super.determinant();

    // Event data
    this->event_data = std::make_shared<NfoldEventData>(
        this->system, state, occ_location, semigrand_canonical_swaps,
        this->potential);

    // Nfold data
    monte::Conversions const &convert =
        get_index_conversions(*this->system, state);
    Index n_allowed_per_unitcell =
        get_n_allowed_per_unitcell(convert, semigrand_canonical_swaps);
    this->nfold_data.n_events_possible =
        static_cast<double>(n_unitcells) * n_allowed_per_unitcell;
  }

  // Make selector
  lotto::RejectionFreeEventSelector event_selector(
      this->event_data->event_calculator,
      clexmonte::make_complete_event_id_list(n_unitcells,
                                             this->event_data->prim_event_list),
      this->event_data->event_list.impact_table,
      std::make_shared<lotto::RandomGenerator>(run_manager.engine));

  // Used to apply selected events: EventID -> monte::OccEvent
  auto get_event_f = [&](EventID const &selected_event_id) {
    // returns a monte::OccEvent
    return this->event_data->event_list.events.at(selected_event_id).event;
  };

  // Run nfold-way
  monte::nfold<EventID>(state, occ_location, this->nfold_data, event_selector,
                        get_event_f, run_manager);
}

}  // namespace nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
