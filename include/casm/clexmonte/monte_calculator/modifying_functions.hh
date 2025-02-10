#ifndef CASM_clexmonte_monte_calculator_modifying_functions
#define CASM_clexmonte_monte_calculator_modifying_functions

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/system/System.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

namespace monte_calculator {

/// \brief Make a state modifying function that sets `mol_composition` and
///     `param_composition` conditions to match the param composition of the
///     state
inline StateModifyingFunction make_match_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return StateModifyingFunction(
      "match.composition",
      "Set `mol_composition` and `param_composition` conditions to match the "
      "calculated composition of the configuration",
      [calculation](state_type &state, monte::OccLocation *occ_location) {
        auto const &system = *calculation->system();

        Eigen::VectorXi const &occupation = get_occupation(state);
        Eigen::VectorXd mol_composition =
            get_composition_calculator(system).mean_num_each_component(
                occupation);
        Eigen::VectorXd param_composition =
            get_composition_converter(system).param_composition(
                mol_composition);
        state.conditions.vector_values["mol_composition"] = mol_composition;
        state.conditions.vector_values["param_composition"] = param_composition;
      });
}

/// \brief Make a state modifying function that enforces the configuration's
///     composition to match the `mol_composition` and/or `param_composition`
///     conditions
inline StateModifyingFunction make_enforce_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return StateModifyingFunction(
      "enforce.composition",
      "Enforce configuration to match `mol_composition` and/or "
      "`param_composition` conditions. Comparison is made using the "
      "\"mol_composition_tol\" calculation parameter, if present.",
      [calculation](state_type &state, monte::OccLocation *occ_location) {
        auto const &system = *calculation->system();
        auto state_data = calculation->state_data();
        if (!state_data) {
          throw std::runtime_error(
              "Error in `enforce.composition`: no state_data");
        }

        if (calculation->params().contains("mol_composition_tol")) {
          if (!calculation->params()["mol_composition_tol"].is_float()) {
            throw std::runtime_error(
                "Error in `enforce.composition`: \"mol_composition_tol\" "
                "parameter is not float");
          }
          calculation->params()["mol_composition_tol"].get<double>();
        }

        ParentInputParser parser{calculation->params()};
        double mol_composition_tol = CASM::TOL;
        parser.optional(mol_composition_tol, "mol_composition_tol");

        std::stringstream ss;
        ss << "Error in `enforce.composition`: error reading calculation "
              "parameters.";
        std::runtime_error error_if_invalid{ss.str()};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

        // Need to check for an OccLocation
        std::unique_ptr<monte::OccLocation> tmp;
        make_temporary_if_necessary(state, occ_location, tmp, *calculation);

        /// - If both present and not consistent, set param_composition to be
        ///   consistent with mol_composition and print warning
        /// - If only one set, set the other to be consistent
        enforce_composition_consistency(
            state, get_composition_converter(system), mol_composition_tol);

        Eigen::VectorXd target_mol_composition =
            get_mol_composition(system, state.conditions);
        monte::RandomNumberGenerator<BaseMonteCalculator::engine_type>
            random_number_generator(calculation->engine());
        clexmonte::enforce_composition(get_occupation(state),
                                       target_mol_composition,
                                       get_composition_calculator(system),
                                       get_semigrand_canonical_swaps(system),
                                       *occ_location, random_number_generator);
      });
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
