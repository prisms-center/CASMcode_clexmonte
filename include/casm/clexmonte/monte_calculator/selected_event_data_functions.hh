#ifndef CASM_clexmonte_monte_calculator_selected_event_data_functions
#define CASM_clexmonte_monte_calculator_selected_event_data_functions

#include "casm/clexmonte/definitions.hh"
#include "casm/monte/sampling/SelectedEventData.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

namespace monte_calculator {

/// \brief Make selected event count collecting function
/// ("selected_event.by_type")
monte::DiscreteVectorIntHistogramFunction make_selected_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event count collecting function
/// ("selected_event.by_equivalent_index")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event count collecting function
/// ("selected_event.by_prim_event_index")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_prim_event_index_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event count collecting functions
/// ("selected_event.<event_type>.by_equivalent_index")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_selected_event_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
