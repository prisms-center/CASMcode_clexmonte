#ifndef CASM_clexmonte_state_modifying_functions
#define CASM_clexmonte_state_modifying_functions

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/composition/CompositionCalculator.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {

// ---
// These methods are used to construct sampling functions. They are templated
// so that they can be reused. The definition documentation should
// state interface requirements for the methods to be applicable and usable in
// a particular context.
//
// Example requirements are:
// - that a conditions `monte::ValueMap` contains scalar "temperature"
// - that the method `ClexData &get_clex(SystemType &,
//   StateType const &, std::string const &key)`
//   exists for template type `SystemType` (i.e. when
//   SystemType=clexmonte::System).
// ---

/// \brief Make a state modifying function that sets `mol_composition`
///     condition equal to the mol composition of the state
///
/// Requires :
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename CalculationType>
StateModifyingFunction make_set_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Make a state modifying function that sets `mol_composition`
///     condition equal to the mol composition of the state
///
/// Requires :
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename CalculationType>
StateModifyingFunction make_set_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return StateModifyingFunction(
      "set_mol_composition",
      "Set `mol_composition` conditions equal to the mol composition of the "
      "state",
      [calculation](state_type &state, monte::OccLocation *occ_location) {
        Eigen::VectorXi const &occupation = get_occupation(state);
        state.conditions.vector_values["mol_composition"] =
            get_composition_calculator(*calculation->system)
                .mean_num_each_component(occupation);
      });
}

}  // namespace clexmonte
}  // namespace CASM

#endif
