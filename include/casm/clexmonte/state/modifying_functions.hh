#ifndef CASM_clexmonte_state_modifying_functions
#define CASM_clexmonte_state_modifying_functions

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/monte/run_management/StateModifyingFunction.hh"

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
state_modifying_function_type make_set_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Make potential energy sampling function ("potential_energy")
///
/// Notes:
/// - This version reads from state.properties, so it works
///   for methods such as `Canonical` and `SemiGrandCanonical`
///   which keep the potential_energy updated, but not `Kinetic`.
///   For `Kinetic` use `make_canonical_potential_energy_f`.
///
/// Requires:
/// - "potential_energy" is a scalar state property

/// \brief Make a state modifying function that sets `mol_composition`
///     condition equal to the mol composition of the state
///
/// Requires :
/// - `composition::CompositionCalculator const &
///   get_composition_calculator(SystemType &)`
template <typename CalculationType>
state_modifying_function_type make_set_mol_composition_f(
    std::shared_ptr<CalculationType> const &calculation) {
  return state_modifying_function_type(
      "set_mol_composition",
      "Set `mol_composition` conditions equal to the mol composition of the "
      "state",
      [calculation](state_type &state) {
        Eigen::VectorXi const &occupation = get_occupation(state);
        state.conditions.vector_values["mol_composition"] =
            get_composition_calculator(*calculation->system)
                .mean_num_each_component(occupation);
      });
}

}  // namespace clexmonte
}  // namespace CASM

#endif
