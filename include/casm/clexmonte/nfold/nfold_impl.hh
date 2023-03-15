#ifndef CASM_clexmonte_nfold_impl
#define CASM_clexmonte_nfold_impl

#include "casm/clexmonte/nfold/nfold.hh"
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/methods/nfold.hh"

namespace CASM {
namespace clexmonte {
namespace nfold {

template <typename EngineType>
Nfold<EngineType>::Nfold(std::shared_ptr<system_type> _system,
                         std::shared_ptr<EngineType> _random_number_engine)
    : semi_grand_canonical::SemiGrandCanonical<EngineType>(
          _system, _random_number_engine) {}

/// \brief Perform a single run, evolving current state
///
/// Notes:
/// - state and occ_location are evolved and end in modified states
template <typename EngineType>
void Nfold<EngineType>::run(state_type &state, monte::OccLocation &occ_location,
                            run_manager_type &run_manager) {
  std::cout << "Nfold::run begin" << std::endl;

  if (!state.conditions.scalar_values.count("temperature")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `temperature` not set.");
  }
  if (!state.conditions.vector_values.count("param_chem_pot")) {
    throw std::runtime_error(
        "Error in Canonical::run: state `param_chem_pot` conditions not set.");
  }

  // Store state info / pointers
  std::cout << "Nfold::run 1" << std::endl;
  this->state = &state;
  this->occ_location = &occ_location;
  this->conditions = make_conditions(*this->system, state);
  Index n_unitcells = this->transformation_matrix_to_super.determinant();

  // Construct potential
  std::cout << "Nfold::run 2" << std::endl;
  typedef semi_grand_canonical::SemiGrandCanonicalPotential potential_type;
  auto potential = std::make_shared<potential_type>(this->system);
  potential->set(this->state, this->conditions);

  // Construct swaps
  std::cout << "Nfold::run 3" << std::endl;
  monte::Conversions const &convert =
      get_index_conversions(*this->system, state);
  monte::OccCandidateList const &occ_candidate_list =
      get_occ_candidate_list(*this->system, state);
  std::vector<monte::OccSwap> grand_canonical_swaps =
      make_grand_canonical_swaps(convert, occ_candidate_list);

  // if same supercell
  // -> just re-set potential & avoid re-constructing event list
  std::cout << "Nfold::run 4" << std::endl;
  if (this->transformation_matrix_to_super ==
          get_transformation_matrix_to_super(state) &&
      this->conditions != nullptr) {
    std::cout << "Nfold::run 5" << std::endl;
    this->event_data->event_calculator->potential = potential;
  } else {
    std::cout << "Nfold::run 6" << std::endl;
    this->transformation_matrix_to_super =
        get_transformation_matrix_to_super(state);
    n_unitcells = this->transformation_matrix_to_super.determinant();

    // Event data
    std::cout << "Nfold::run 7" << std::endl;
    this->event_data = std::make_shared<NfoldEventData>(
        this->system, state, occ_location, grand_canonical_swaps, potential);

    // Nfold data
    std::cout << "Nfold::run 8" << std::endl;
    Index n_allowed_per_unitcell =
        get_n_allowed_per_unitcell(convert, grand_canonical_swaps);
    this->nfold_data.n_events_possible =
        static_cast<double>(n_unitcells) * n_allowed_per_unitcell;
  }

  // Make selector
  std::cout << "Nfold::run 9" << std::endl;
  lotto::RejectionFreeEventSelector event_selector(
      this->event_data->event_calculator,
      clexmonte::make_complete_event_id_list(n_unitcells,
                                             this->event_data->prim_event_list),
      this->event_data->event_list.impact_table);

  // Used to apply selected events: EventID -> monte::OccEvent
  std::cout << "Nfold::run 10" << std::endl;
  auto get_event_f = [&](EventID const &selected_event_id) {
    // returns a monte::OccEvent
    return this->event_data->event_list.events.at(selected_event_id).event;
  };

  // Run nfold-way
  std::cout << "Nfold::run 11" << std::endl;
  monte::nfold(state, occ_location, this->nfold_data, event_selector,
               get_event_f, run_manager);
  std::cout << "Nfold::run end" << std::endl;
}

}  // namespace nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
