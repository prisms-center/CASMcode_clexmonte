#ifndef CASM_clexmonte_monte_calculator_kinetic_sampling_functions
#define CASM_clexmonte_monte_calculator_kinetic_sampling_functions

#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/misc/diffusion_calculations.hh"
#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/Correlations.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

namespace CASM {
namespace clexmonte {

class MonteCalculator;

namespace monte_calculator {

/// \brief Make center of mass isotropic squared displacement sampling function
///     ("mean_R_squared_collective_isotropic")
state_sampling_function_type make_mean_R_squared_collective_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make center of mass anisotropic squared displacement sampling
/// function
/// ("mean_R_squared_collective_anisotropic")
state_sampling_function_type make_mean_R_squared_collective_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make tracer isotropic squared displacement sampling function
///     ("mean_R_squared_individual_isotropic")
state_sampling_function_type make_mean_R_squared_individual_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make tracer anisotropic squared displacement sampling function
///     ("mean_R_squared_individual_anisotropic")
state_sampling_function_type make_mean_R_squared_individual_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make isotropic Onsager kinetic coefficient sampling function
///     ("L_isotropic")
state_sampling_function_type make_L_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make anisotropic Onsager kinetic coefficient sampling function
///     ("L_anisotropic")
state_sampling_function_type make_L_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make isotropic tracer diffusion coefficient sampling function
///     ("D_tracer_isotropic")
state_sampling_function_type make_D_tracer_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make anisotropic tracer diffusion coefficient sampling function
///     ("D_tracer_anisotropic")
state_sampling_function_type make_D_tracer_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make delta_n_jumps(i) / n_atoms(i) ("jumps_per_atom_by_type")
state_sampling_function_type make_jumps_per_atom_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make delta_n_jumps(i) / delta_n_events ("jumps_per_event_by_type")
state_sampling_function_type make_jumps_per_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make delta_n_jumps(i) / n_atoms(i) / delta_n_events
/// ("jumps_per_atom_per_event_by_type")
state_sampling_function_type make_jumps_per_atom_per_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

// --- Histogram sampling functions ---

/// \brief Sample the change in histogram counts or fraction of the total
/// change in histogram counts during the last sampling period
template <typename ValueType, typename CompareType, typename HistogramType>
struct HistogramSamplingFunctionT {
  /// The calculation pointer
  std::shared_ptr<MonteCalculator> calculation;

  /// The function name
  std::string sampling_function_name;

  /// The histogram name
  std::string histogram_name;

  /// If true, sample the delta count of the histogram, otherwise sample the
  /// fraction
  bool sample_count;

  /// The values to sample -> {label, index}, where
  /// - label: corresponds to sampled vector output component name
  /// - index: the sampled vector component index
  std::map<ValueType, std::pair<std::string, Index>, CompareType> values;

  /// The previous count of the specified values
  Eigen::VectorXd prev_count;

  /// \brief Constructor
  HistogramSamplingFunctionT(
      std::shared_ptr<MonteCalculator> const &_calculation,
      std::string const &_sampling_function_name,
      std::string const &_histogram_name, bool _sample_count,
      std::map<ValueType, std::string, CompareType> _value_labels);

  /// The number of values sampled (+ 1 for "other" values)
  Index size() const { return values.size() + 1; }

  /// The shape of the output (={size()})
  std::vector<Index> shape() const { return {size()}; }

  /// The component names (labels, in order determined by value_labels)
  std::vector<std::string> component_names() const;

  Eigen::VectorXd operator()();
};

typedef HistogramSamplingFunctionT<Eigen::VectorXl,
                                   monte::LexicographicalCompare,
                                   monte::DiscreteVectorIntHistogram>
    DiscreteVectorIntHistogramSamplingFunction;

typedef HistogramSamplingFunctionT<Eigen::VectorXd,
                                   monte::FloatLexicographicalCompare,
                                   monte::DiscreteVectorFloatHistogram>
    DiscreteVectorFloatHistogramSamplingFunction;

/// Get a histogram of a particular type from calculation->selected_event_data
template <typename HistogramType>
HistogramType get_histogram(std::shared_ptr<MonteCalculator> const &calculation,
                            std::string sampling_function_name,
                            std::string histogram_name);

/// Get a DiscreteVectorIntHistogram from calculation->selected_event_data
template <>
monte::DiscreteVectorIntHistogram
get_histogram<monte::DiscreteVectorIntHistogram>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name);

/// Get a DiscreteVectorFloatHistogram from calculation->selected_event_data
template <>
monte::DiscreteVectorFloatHistogram
get_histogram<monte::DiscreteVectorFloatHistogram>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name);

/// Get a PartitionedHistogram1D from calculation->selected_event_data
template <>
monte::PartitionedHistogram1D get_histogram<monte::PartitionedHistogram1D>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name);

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_type")
state_sampling_function_type make_selected_event_count_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_type")
state_sampling_function_type make_selected_event_fraction_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_equivalent_index")
state_sampling_function_type make_selected_event_count_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_equivalent_index")
state_sampling_function_type make_selected_event_fraction_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_equivalent_index_and_direction")
state_sampling_function_type
make_selected_event_count_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_equivalent_index_and_direction")
state_sampling_function_type
make_selected_event_fraction_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.count.<event_type>.by_equivalent_index")
std::vector<state_sampling_function_type>
make_selected_event_count_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.<event_type>.by_equivalent_index")
std::vector<state_sampling_function_type>
make_selected_event_fraction_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation);

// -- Inline definitions --

/// Get a histogram of a particular type from calculation->selected_event_data
template <typename HistogramType>
HistogramType get_histogram(std::shared_ptr<MonteCalculator> const &calculation,
                            std::string sampling_function_name,
                            std::string histogram_name) {
  throw std::runtime_error(
      "Error in " + sampling_function_name +
      " sampling function: not of a supported histogram type");
}

/// \brief Constructor
template <typename ValueType, typename CompareType, typename HistogramType>
HistogramSamplingFunctionT<ValueType, CompareType, HistogramType>::
    HistogramSamplingFunctionT(
        std::shared_ptr<MonteCalculator> const &_calculation,
        std::string const &_sampling_function_name,
        std::string const &_histogram_name, bool _sample_count,
        std::map<ValueType, std::string, CompareType> _value_labels)
    : calculation(_calculation),
      sampling_function_name(_sampling_function_name),
      histogram_name(_histogram_name),
      sample_count(_sample_count),
      values(_value_labels.key_comp()) {
  Index i = 0;
  for (auto const &pair : _value_labels) {
    values.emplace(pair.first, std::make_pair(pair.second, i));
    ++i;
  }

  // Set `prev_count` to 0s
  prev_count = Eigen::VectorXd::Zero(size());
}

/// The component names (labels, in order determined by value_labels)
template <typename ValueType, typename CompareType, typename HistogramType>
std::vector<std::string> HistogramSamplingFunctionT<
    ValueType, CompareType, HistogramType>::component_names() const {
  std::vector<std::string> _component_names;
  for (auto const &pair : values) {
    _component_names.push_back(pair.second.first);
  }
  _component_names.push_back("other");
  return _component_names;
}

/// \brief Call operator
template <typename ValueType, typename CompareType, typename HistogramType>
Eigen::VectorXd HistogramSamplingFunctionT<ValueType, CompareType,
                                           HistogramType>::operator()() {
  HistogramType const &hist = get_histogram<HistogramType>(
      calculation, sampling_function_name, histogram_name);

  // Set current count to 0s
  Eigen::VectorXd current_count = Eigen::VectorXd::Zero(size());

  // Add counts for each value
  for (auto const &pair : hist.value_counts()) {
    auto it = values.find(pair.first);
    if (it != values.end()) {
      current_count(it->second.second) += pair.second;
    } else {
      // "other" values lumped together
      current_count(size() - 1) += pair.second;
    }
  }

  // "other" values lumped together includes out-of-range values
  current_count(size() - 1) += hist.out_of_range_count();

  // Calculate delta count
  Eigen::VectorXd delta_count = current_count - prev_count;

  // Set prev_count to current_count
  prev_count = current_count;

  // Return delta count or fraction
  if (sample_count) {
    return delta_count;
  } else {
    return delta_count / delta_count.sum();
  }
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM

#endif
