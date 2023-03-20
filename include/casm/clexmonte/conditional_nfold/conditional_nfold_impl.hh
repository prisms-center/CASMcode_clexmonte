#ifndef CASM_clexmonte_conditional_nfold_impl
#define CASM_clexmonte_conditional_nfold_impl

#include "casm/clexmonte/conditional_nfold/conditional_nfold.hh"
#include "casm/clexmonte/nfold/nfold_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/nfold.hh"

namespace CASM {
namespace clexmonte {
namespace conditional_nfold {

template <typename EngineType>
ConditionalNfold<EngineType>::ConditionalNfold(
    std::shared_ptr<system_type> _system,
    std::shared_ptr<EngineType> _random_number_engine,
    double _nfold_begin_acceptance_rate, double _nfold_end_acceptance_rate)
    : nfold::Nfold<EngineType>(_system, _random_number_engine),
      nfold_begin_acceptance_rate(_nfold_begin_acceptance_rate),
      nfold_end_acceptance_rate(_nfold_end_acceptance_rate) {}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void ConditionalNfold<EngineType>::run(
    state_type &state, monte::OccLocation &occ_location,
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
  this->transformation_matrix_to_super =
      get_transformation_matrix_to_super(state);
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);
  Index n_unitcells = this->transformation_matrix_to_super.determinant();

  // Construct potential
  typedef semi_grand_canonical::SemiGrandCanonicalPotential potential_type;
  auto potential = std::make_shared<potential_type>(this->system);
  potential->set(this->state, this->conditions);

  // Construct swaps
  monte::Conversions const &convert =
      get_index_conversions(*this->system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*this->system, state);
  std::vector<monte::OccSwap> grand_canonical_swaps =
      make_grand_canonical_swaps(convert, occ_candidate_list);

  // Event data
  this->event_data = std::make_shared<nfold::NfoldEventData>(
      this->system, state, occ_location, grand_canonical_swaps, potential);

  // Nfold data
  Index n_allowed_per_unitcell =
      get_n_allowed_per_unitcell(convert, grand_canonical_swaps);
  this->nfold_data.n_events_possible =
      static_cast<double>(n_unitcells) * n_allowed_per_unitcell;

  // Make selector
  lotto::RejectionFreeEventSelector event_selector(
      this->event_data->event_calculator,
      clexmonte::make_complete_event_id_list(n_unitcells,
                                             this->event_data->prim_event_list),
      this->event_data->event_list.impact_table);

  // Used to apply selected events: EventID -> monte::OccEvent
  auto get_event_f = [&](EventID const &selected_event_id) {
    // returns a monte::OccEvent
    return this->event_data->event_list.events.at(selected_event_id).event;
  };

  // Run nfold-way
  monte::nfold(state, occ_location, this->nfold_data, event_selector,
               get_event_f, run_manager);
}

}  // namespace conditional_nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
