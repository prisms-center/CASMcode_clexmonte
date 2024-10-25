#ifndef CASM_clexmonte_monte_calculator_selected_event_functions
#define CASM_clexmonte_monte_calculator_selected_event_functions

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/LocalOrbitCompositionCalculator.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

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
/// ("selected_event.by_equivalent_index_and_direction")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event count collecting functions
/// ("selected_event.<event_type>.by_equivalent_index")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_selected_event_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

// -- Local orbit composition --

/// \brief Helper class for collecting local orbit composition data
///
/// Generally, the particular state being calculated is not known until the
/// MonteCalculator `run` method is called. This class manages checking the
/// MonteCalculator state and constructing a LocalOrbitCompositionCalculator
/// for that state.
///
/// Notes:
/// - Holds a MonteCalculator pointer and LocalOrbitCompositionCalculatorData
///   which is used to construct a LocalOrbitCompositionCalculator for the
///   state currently being calculated.
/// - After the LocalOrbitCompositionCalculator is constructed, the `value`
///   method is used to calculate the local orbit composition for a
///   particular event in the state.
/// - Assumes that the MonteCalculator state's supercell does not change during
///   the lifetime of the LocalOrbitCompositionCollector.
class LocalOrbitCompositionCollector {
 public:
  LocalOrbitCompositionCollector(
      std::shared_ptr<MonteCalculator> calculation,
      std::shared_ptr<LocalOrbitCompositionCalculatorData const> data);

  /// \brief Return the shape of the output values
  std::vector<Index> shape() const;

  /// \brief Return names for value components (col-major unrolling)
  std::vector<std::string> component_names() const;

  /// \brief Evaluate the local orbit composition for a particular event
  Eigen::MatrixXi const &value(Index unitcell_index, Index equivalent_index);

  /// \brief Get the MonteCalculator being used
  std::shared_ptr<MonteCalculator> calculation() const { return m_calculation; }

  /// \brief Get the LocalOrbitCompositionCalculatorData being used
  std::shared_ptr<LocalOrbitCompositionCalculatorData const> data() const {
    return m_data;
  }

  /// \brief Get the LocalOrbitCompositionCalculator being used
  std::shared_ptr<LocalOrbitCompositionCalculator>
  local_orbit_composition_calculator() const {
    return m_local_orbit_composition_calculator;
  }

  /// \brief Reset the LocalOrbitCompositionCalculator using the current state
  /// being calculated by the MonteCalculator
  void reset();

 private:
  std::shared_ptr<MonteCalculator> m_calculation;
  std::shared_ptr<LocalOrbitCompositionCalculatorData const> m_data;
  std::shared_ptr<LocalOrbitCompositionCalculator>
      m_local_orbit_composition_calculator;
};

/// \brief Make local orbit composition collecting functions
/// ("local_orbit_composition.<key>")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_local_orbit_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation);

// -- dE_activated --

/// \brief Make dE_activated collecting function
/// ("dE_activated.by_type")
monte::PartitionedHistogramFunction<double> make_dE_activated_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make dE_activated collecting function
/// ("dE_activated.by_equivalent_index")
monte::PartitionedHistogramFunction<double>
make_dE_activated_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation);

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
