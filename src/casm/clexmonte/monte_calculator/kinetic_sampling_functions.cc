#include "casm/clexmonte/monte_calculator/kinetic_sampling_functions.hh"

#include "casm/clexmonte/canonical/canonical.hh"
#include "casm/clexmonte/misc/diffusion_calculations.hh"
#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/Correlations.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/sampling/SelectedEventFunctions.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {
namespace monte_calculator {

namespace {

// --- Helper functions ---

// Should be OK to evaluate after a calculator has had
// selected_event_functions added
monte::DiscreteVectorIntHistogramFunction _get_vector_int_histogram_function(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &functions =
      calculation->selected_event_functions()->discrete_vector_int_functions;
  auto it = functions.find(histogram_name);
  if (it == functions.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "Selected event function '" +
                             histogram_name + "' not found");
  }
  return it->second;
}

// Should be OK to evaluate after a calculator has had
// selected_event_functions added
monte::DiscreteVectorFloatHistogramFunction
_get_vector_float_histogram_function(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &functions =
      calculation->selected_event_functions()->discrete_vector_float_functions;
  auto it = functions.find(histogram_name);
  if (it == functions.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "Selected event function '" +
                             histogram_name + "' not found");
  }
  return it->second;
}

// Should be OK to evaluate after a calculator has had
// selected_event_functions added
monte::PartitionedHistogramFunction<double>
_get_continuous_1d_histogram_function(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &functions =
      calculation->selected_event_functions()->continuous_1d_functions;
  auto it = functions.find(histogram_name);
  if (it == functions.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "Selected event function '" +
                             histogram_name + "' not found");
  }
  return it->second;
}

// Should be OK to evaluate after a run begins,
// when evaluating a sampling function
monte::SelectedEventData const &_get_selected_event_data(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name) {
  if (!calculation->selected_event_data()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "selected_event_data is null");
  }
  return *calculation->selected_event_data();
}

// Should be OK to evaluate after a run begins
monte::DiscreteVectorIntHistogram const &_get_vector_int_histogram(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &selected_event_data =
      _get_selected_event_data(calculation, sampling_function_name);
  auto const &histograms = selected_event_data.discrete_vector_int_histograms;
  auto it = histograms.find(histogram_name);
  if (it == histograms.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "selected event data '" +
                             histogram_name + "' is not being collected");
  }
  return it->second;
}

// Should be OK to evaluate after a run begins
monte::DiscreteVectorFloatHistogram const &_get_vector_float_histogram(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &selected_event_data =
      _get_selected_event_data(calculation, sampling_function_name);
  auto const &histograms = selected_event_data.discrete_vector_float_histograms;
  auto it = histograms.find(histogram_name);
  if (it == histograms.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "selected event data '" +
                             histogram_name + "' is not being collected");
  }
  return it->second;
}

// Should be OK to evaluate after a run begins
monte::PartitionedHistogram1D const &_get_continuous_1d_histogram(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  auto const &selected_event_data =
      _get_selected_event_data(calculation, sampling_function_name);
  auto const &histograms = selected_event_data.continuous_1d_histograms;
  auto it = histograms.find(histogram_name);
  if (it == histograms.end()) {
    throw std::runtime_error("Error in " + sampling_function_name +
                             " sampling function: "
                             "selected event data '" +
                             histogram_name + "' is not being collected");
  }
  return it->second;
}

}  // namespace

/// Get a DiscreteVectorIntHistogram from calculation->selected_event_data
template <>
monte::DiscreteVectorIntHistogram
get_histogram<monte::DiscreteVectorIntHistogram>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  return _get_vector_int_histogram(calculation, sampling_function_name,
                                   histogram_name);
}

/// Get a DiscreteVectorFloatHistogram from calculation->selected_event_data
template <>
monte::DiscreteVectorFloatHistogram
get_histogram<monte::DiscreteVectorFloatHistogram>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  return _get_vector_float_histogram(calculation, sampling_function_name,
                                     histogram_name);
}

/// Get a PartitionedHistogram1D from calculation->selected_event_data
template <>
monte::PartitionedHistogram1D get_histogram<monte::PartitionedHistogram1D>(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string sampling_function_name, std::string histogram_name) {
  return _get_continuous_1d_histogram(calculation, sampling_function_name,
                                      histogram_name);
}

/// \brief Make center of mass isotropic squared displacement sampling function
///     ("mean_R_squared_collective_isotropic")
state_sampling_function_type make_mean_R_squared_collective_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<CollectiveIsotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "mean_R_squared_collective_isotropic",
      R"(Samples \frac{1}{N} \left(\sum_\zeta \Delta R^\zeta_{i} \right) \dot \left(\sum_\zeta \Delta R^\zeta_{j} \right))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        Eigen::VectorXd result = mean_R_squared_collective_isotropic(
            name_list, name_index_list, delta_R);
        return result;
      });
}

/// \brief Make center of mass anisotropic squared displacement sampling
/// function
/// ("mean_R_squared_collective_anisotropic")
state_sampling_function_type make_mean_R_squared_collective_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<CollectiveAnisotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "mean_R_squared_collective_anisotropic",
      R"(Samples \frac{1}{N} \left(\sum_\zeta \Delta R^\zeta_{i,\alpha} \right) \left(\sum_\zeta \Delta R^\zeta_{j,\beta} \right))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        return mean_R_squared_collective_anisotropic(name_list, name_index_list,
                                                     delta_R);
      });
}

/// \brief Make tracer isotropic squared displacement sampling function
///     ("mean_R_squared_individual_isotropic")
state_sampling_function_type make_mean_R_squared_individual_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<IndividualIsotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "mean_R_squared_individual_isotropic",
      R"(Samples \frac{1}{N_i} \sum_\zeta \left(\Delta R^\zeta_{i} \dot \Delta R^\zeta_{i}\right))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        return mean_R_squared_individual_isotropic(name_list, name_index_list,
                                                   delta_R);
      });
}

/// \brief Make tracer anisotropic squared displacement sampling function
///     ("mean_R_squared_individual_anisotropic")
state_sampling_function_type make_mean_R_squared_individual_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<IndividualAnisotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "mean_R_squared_individual_anisotropic",  // individual
      R"(Samples \frac{1}{N_i} \sum_\zeta \left(\Delta R^\zeta_{i,\alpha} \Delta R^\zeta_{i,\beta}\right))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        return mean_R_squared_individual_anisotropic(name_list, name_index_list,
                                                     delta_R);
      });
}

/// \brief Make isotropic Onsager kinetic coefficient sampling function
///     ("L_isotropic")
state_sampling_function_type make_L_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<CollectiveIsotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "L_isotropic",
      R"(Samples \frac{1}{N} \left(\sum_\zeta \Delta R^\zeta_{i} \right) \dot \left(\sum_\zeta \Delta R^\zeta_{j} \right) / (2 d \Delta t))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        auto const &time_curr = kmc_data.time;
        auto const &time_prev =
            kmc_data.prev_time.at(kmc_data.sampling_fixture_label);
        double delta_time = time_curr - time_prev;

        double dim = system.n_dimensions;
        double normalization = (2.0 * dim * delta_time);

        Eigen::VectorXd mean_R_squared = mean_R_squared_collective_isotropic(
            name_list, name_index_list, delta_R);
        Eigen::VectorXd L_anisotropic = mean_R_squared / normalization;
        return L_anisotropic;
      });
}

/// \brief Make anisotropic Onsager kinetic coefficient sampling function
///     ("L_anisotropic")
state_sampling_function_type make_L_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<CollectiveAnisotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "L_anisotropic",
      R"(Samples \frac{1}{N} \left(\sum_\zeta \Delta R^\zeta_{i} \right) \dot \left(\sum_\zeta \Delta R^\zeta_{j} \right) / (2 \Delta t))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        auto const &time_curr = kmc_data.time;
        auto const &time_prev =
            kmc_data.prev_time.at(kmc_data.sampling_fixture_label);
        double delta_time = time_curr - time_prev;

        double normalization = (2.0 * delta_time);

        Eigen::VectorXd mean_R_squared = mean_R_squared_collective_anisotropic(
            name_list, name_index_list, delta_R);
        Eigen::VectorXd L_anisotropic = mean_R_squared / normalization;
        return L_anisotropic;
      });
}

/// \brief Make isotropic tracer diffusion coefficient sampling function
///     ("D_tracer_isotropic")
state_sampling_function_type make_D_tracer_isotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<IndividualIsotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "D_tracer_isotropic",
      R"(Samples \frac{1}{N_i} \sum_\zeta \left(\Delta R^\zeta_{i} \dot \Delta R^\zeta_{i}\right) / (2 d \Delta t))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        auto const &time_curr = kmc_data.time;
        auto const &time_prev =
            kmc_data.prev_time.at(kmc_data.sampling_fixture_label);
        double delta_time = time_curr - time_prev;

        double dim = system.n_dimensions;
        double normalization = (2.0 * dim * delta_time);

        Eigen::VectorXd mean_R_squared = mean_R_squared_individual_isotropic(
            name_list, name_index_list, delta_R);
        Eigen::VectorXd D_tracer_isotropic = mean_R_squared / normalization;
        return D_tracer_isotropic;
      });
}

/// \brief Make anisotropic tracer diffusion coefficient sampling function
///     ("D_tracer_anisotropic")
state_sampling_function_type make_D_tracer_anisotropic_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  auto const &name_list = event_system->atom_name_list;
  std::vector<std::string> component_names =
      make_component_names<IndividualAnisotropicCounter>(name_list);

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "D_tracer_anisotropic",
      R"(Samples \frac{1}{N_i} \sum_\zeta \left(\Delta R^\zeta_{i} \dot \Delta R^\zeta_{i}\right) / (2 \Delta t))",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = *calculation->system();
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = *calculation->kmc_data();
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &R_curr = kmc_data.atom_positions_cart;
        auto const &R_prev = kmc_data.prev_atom_positions_cart.at(
            kmc_data.sampling_fixture_label);
        Eigen::MatrixXd delta_R = R_curr - R_prev;

        auto const &time_curr = kmc_data.time;
        auto const &time_prev =
            kmc_data.prev_time.at(kmc_data.sampling_fixture_label);
        double delta_time = time_curr - time_prev;

        double normalization = (2.0 * delta_time);

        Eigen::VectorXd mean_R_squared = mean_R_squared_individual_anisotropic(
            name_list, name_index_list, delta_R);
        Eigen::VectorXd D_tracer_anisotropic = mean_R_squared / normalization;
        return D_tracer_anisotropic;
      });
}

/// \brief Make delta_n_jumps(i) / n_atoms(i) ("jumps_per_atom_by_type")
state_sampling_function_type make_jumps_per_atom_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  std::vector<std::string> component_names = event_system->atom_name_list;

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  std::shared_ptr<Index> prev_n_events = std::make_shared<Index>(0);
  std::shared_ptr<Eigen::VectorXd> prev_sum_n_jumps =
      std::make_shared<Eigen::VectorXd>(
          Eigen::VectorXd::Zero(component_names.size()));

  return state_sampling_function_type(
      "jumps_per_atom_by_type",  // individual
      R"(Mean number of jumps per atom for each atom type over the last sampling period)",
      component_names,  // component names
      shape, [calculation, prev_n_events, prev_sum_n_jumps]() {
        auto const &system = *calculation->system();
        auto const &kmc_data = *calculation->kmc_data();
        auto const &occ_location = *calculation->state_data()->occ_location;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &n_jumps = occ_location.current_atom_n_jumps();

        auto const &sampling_fixture = *kmc_data.sampling_fixture;
        auto const &counter = sampling_fixture.counter();
        double steps_per_pass = counter.steps_per_pass;
        double step = counter.step;
        double pass = counter.pass;
        double n_events = steps_per_pass * pass + step;

        // reset stored data if necessary
        if (*prev_n_events > n_events) {
          *prev_sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
          *prev_n_events = 0;
        }

        Eigen::VectorXd n_atoms = Eigen::VectorXd::Zero(name_list.size());
        Eigen::VectorXd sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
        for (Index i = 0; i < n_jumps.size(); ++i) {
          n_atoms(name_index_list[i]) += 1.0;
          sum_n_jumps(name_index_list[i]) += n_jumps[i];
        }
        Eigen::VectorXd delta_n_jumps = sum_n_jumps - *prev_sum_n_jumps;

        // jumps_per_atom_by_type = delta_n_jumps(i) / n_atoms(i);
        Eigen::VectorXd jumps_per_atom_by_type =
            Eigen::VectorXd::Zero(name_list.size());
        for (Index i = 0; i < name_list.size(); ++i) {
          jumps_per_atom_by_type(i) = delta_n_jumps(i) / n_atoms(i);
        }

        *prev_sum_n_jumps = sum_n_jumps;
        *prev_n_events = n_events;

        return jumps_per_atom_by_type;
      });
}

/// \brief Make delta_n_jumps(i) / delta_n_events ("jumps_per_event_by_type")
state_sampling_function_type make_jumps_per_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  std::vector<std::string> component_names = event_system->atom_name_list;

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  std::shared_ptr<Index> prev_n_events = std::make_shared<Index>(0);
  std::shared_ptr<Eigen::VectorXd> prev_sum_n_jumps =
      std::make_shared<Eigen::VectorXd>(
          Eigen::VectorXd::Zero(component_names.size()));

  return state_sampling_function_type(
      "jumps_per_event_by_type",  // individual
      R"(Mean number of jumps per event for each atom type over the last sampling period)",
      component_names,  // component names
      shape, [calculation, prev_n_events, prev_sum_n_jumps]() {
        auto const &system = *calculation->system();
        auto const &kmc_data = *calculation->kmc_data();
        auto const &occ_location = *calculation->state_data()->occ_location;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &n_jumps = occ_location.current_atom_n_jumps();

        auto const &sampling_fixture = *kmc_data.sampling_fixture;
        auto const &counter = sampling_fixture.counter();
        double steps_per_pass = counter.steps_per_pass;
        double step = counter.step;
        double pass = counter.pass;
        double n_events = steps_per_pass * pass + step;

        // reset stored data if necessary
        if (*prev_n_events > n_events) {
          *prev_sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
          *prev_n_events = 0;
        }
        double delta_n_events = n_events - *prev_n_events;

        Eigen::VectorXd n_atoms = Eigen::VectorXd::Zero(name_list.size());
        Eigen::VectorXd sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
        for (Index i = 0; i < n_jumps.size(); ++i) {
          n_atoms(name_index_list[i]) += 1.0;
          sum_n_jumps(name_index_list[i]) += n_jumps[i];
        }
        Eigen::VectorXd delta_n_jumps = sum_n_jumps - *prev_sum_n_jumps;

        // jumps_per_event_by_type = delta_n_jumps(i) / delta_n_events;
        Eigen::VectorXd jumps_per_event_by_type =
            delta_n_jumps / delta_n_events;

        *prev_sum_n_jumps = sum_n_jumps;
        *prev_n_events = n_events;

        return jumps_per_event_by_type;
      });
}

/// \brief Make delta_n_jumps(i) / n_atoms(i) / delta_n_events
/// ("jumps_per_atom_per_event_by_type")
state_sampling_function_type make_jumps_per_atom_per_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system();
  auto event_system = get_event_system(system);
  std::vector<std::string> component_names = event_system->atom_name_list;

  std::vector<Index> shape;
  shape.push_back(component_names.size());

  std::shared_ptr<Index> prev_n_events = std::make_shared<Index>(0);
  std::shared_ptr<Eigen::VectorXd> prev_sum_n_jumps =
      std::make_shared<Eigen::VectorXd>(
          Eigen::VectorXd::Zero(component_names.size()));

  return state_sampling_function_type(
      "jumps_per_atom_per_event_by_type",  // individual
      R"(Mean number of jumps per event for each atom type over the last sampling period)",
      component_names,  // component names
      shape, [calculation, prev_n_events, prev_sum_n_jumps]() {
        auto const &system = *calculation->system();
        auto const &kmc_data = *calculation->kmc_data();
        auto const &occ_location = *calculation->state_data()->occ_location;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list = kmc_data.atom_name_index_list;
        auto const &n_jumps = occ_location.current_atom_n_jumps();

        auto const &sampling_fixture = *kmc_data.sampling_fixture;
        auto const &counter = sampling_fixture.counter();
        double steps_per_pass = counter.steps_per_pass;
        double step = counter.step;
        double pass = counter.pass;
        double n_events = steps_per_pass * pass + step;

        // reset stored data if necessary
        if (*prev_n_events > n_events) {
          *prev_sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
          *prev_n_events = 0;
        }
        double delta_n_events = n_events - *prev_n_events;

        Eigen::VectorXd n_atoms = Eigen::VectorXd::Zero(name_list.size());
        Eigen::VectorXd sum_n_jumps = Eigen::VectorXd::Zero(name_list.size());
        for (Index i = 0; i < n_jumps.size(); ++i) {
          n_atoms(name_index_list[i]) += 1.0;
          sum_n_jumps(name_index_list[i]) += n_jumps[i];
        }
        Eigen::VectorXd delta_n_jumps = sum_n_jumps - *prev_sum_n_jumps;

        // jumps_per_atom_per_event_by_type = delta_n_jumps(i) / n_atoms(i) /
        // delta_n_events;
        Eigen::VectorXd jumps_per_atom_per_event_by_type =
            Eigen::VectorXd::Zero(name_list.size());
        for (Index i = 0; i < name_list.size(); ++i) {
          jumps_per_atom_per_event_by_type(i) =
              delta_n_jumps(i) / n_atoms(i) / delta_n_events;
        }

        *prev_sum_n_jumps = sum_n_jumps;
        *prev_n_events = n_events;

        return jumps_per_atom_per_event_by_type;
      });
}

// -- Selected event ----------------------------

// -- Selected event by type --

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_type")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_type` has already been added
state_sampling_function_type make_selected_event_count_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "selected_event.count.by_type";
  std::string desc =
      "Selected event count by event type. Requires "
      "selected_event.by_type selected event data is collected.";

  std::string histogram_name = "selected_event.by_type";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = true;

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_type")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_type` has already been added
state_sampling_function_type make_selected_event_fraction_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "selected_event.fraction.by_type";
  std::string desc =
      "Selected event fraction by event type. Requires "
      "selected_event.by_type selected event data is collected.";

  std::string histogram_name = "selected_event.by_type";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = false;  // sample fraction

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

// -- Selected event by equivalent index --

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_equivalent_index")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_equivalent_index` has already been added
state_sampling_function_type make_selected_event_count_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "selected_event.count.by_equivalent_index";
  std::string desc =
      "Selected event count, for all events by equivalent index. In the set of "
      "symmetrically equivalent events, events with the same equivalent index"
      "differ only by a translation. Requires selected_event.by_type selected "
      "event data is collected.";

  std::string histogram_name = "selected_event.by_equivalent_index";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = true;

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_equivalent_index")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_equivalent_index` has already been added
state_sampling_function_type make_selected_event_fraction_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "selected_event.fraction.by_equivalent_index";
  std::string desc =
      "Selected event count, for all events by equivalent index. In the set of "
      "symmetrically equivalent events, events with the same equivalent index"
      "differ only by a translation. Requires selected_event.by_type selected "
      "event data is collected.";

  std::string histogram_name = "selected_event.by_equivalent_index";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = false;  // sample fraction

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

// -- Selected event by prim event index --

/// \brief Make selected event type sampling function
/// ("selected_event.count.by_equivalent_index_and_direction")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_equivalent_index_and_direction` has already been added
state_sampling_function_type
make_selected_event_count_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "selected_event.count.by_equivalent_index_and_direction";
  std::string desc =
      "Selected event count, for all events by prim event index. In the set of "
      "symmetrically equivalent events, events with the same prim event index"
      "differ only by a translation and jump in the same direction. Requires "
      "selected_event.by_equivalent_index_and_direction selected event data is "
      "collected.";

  std::string histogram_name =
      "selected_event.by_equivalent_index_and_direction";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = true;

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.by_equivalent_index_and_direction")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.by_equivalent_index_and_direction` has already been added
state_sampling_function_type
make_selected_event_fraction_by_equivalent_index_and_direction_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name =
      "selected_event.fraction.by_equivalent_index_and_direction";
  std::string desc =
      "Selected event count, for all events by prim event index. In the set of "
      "symmetrically equivalent events, events with the same prim event index"
      "differ only by a translation and jump in the same direction. Requires "
      "selected_event.by_equivalent_index_and_direction selected event data is "
      "collected.";

  std::string histogram_name =
      "selected_event.by_equivalent_index_and_direction";
  auto const &event_f =
      _get_vector_int_histogram_function(calculation, name, histogram_name);
  bool sample_count = false;  // sample fraction

  DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
      calculation, name, histogram_name, sample_count, *event_f.value_labels);

  return state_sampling_function_type(name, desc,
                                      hist_sampling_f.component_names(),
                                      hist_sampling_f.shape(), hist_sampling_f);
}

// -- Selected event by equivalent index, per event type --

/// \brief Make selected event type sampling function
/// ("selected_event.count.<event_type>.by_equivalent_index")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.<event_type>.by_equivalent_index` has already been added
/// - This is the same as i.e. `selected_event.count.by_equivalent_index` but
///   for a single event type instead of all event types
std::vector<state_sampling_function_type>
make_selected_event_count_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // get names in alphabetical order
  std::set<std::string> keys;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    keys.insert(x.event_type_name);
  }

  std::vector<state_sampling_function_type> f_list;

  for (std::string event_type_name : keys) {
    std::string name =
        "selected_event.count." + event_type_name + ".by_equivalent_index";
    std::string desc =
        "Selected event count, for all events of a single event type by "
        "equivalent index. In the set of symmetrically equivalent events, "
        "events with the same equivalent index differ only by a translation. "
        "Requires selected_event." +
        event_type_name +
        ".by_equivalent_index selected event data is collected.";

    std::string histogram_name =
        "selected_event." + event_type_name + ".by_equivalent_index";
    auto const &event_f =
        _get_vector_int_histogram_function(calculation, name, histogram_name);
    bool sample_count = true;

    DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
        calculation, name, histogram_name, sample_count, *event_f.value_labels);

    f_list.emplace_back(name, desc, hist_sampling_f.component_names(),
                        hist_sampling_f.shape(), hist_sampling_f);
  }
  return f_list;
}

/// \brief Make selected event type sampling function
/// ("selected_event.fraction.<event_type>.by_equivalent_index")
///
/// Notes:
/// - This requires that the Selected event function
///   `selected_event.<event_type>.by_equivalent_index` has already been added
/// - This is the same as i.e. `selected_event.fraction.by_equivalent_index` but
///   for a single event type instead of all event types
std::vector<state_sampling_function_type>
make_selected_event_fraction_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // get names in alphabetical order
  std::set<std::string> keys;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    keys.insert(x.event_type_name);
  }

  std::vector<state_sampling_function_type> f_list;

  for (std::string event_type_name : keys) {
    std::string name =
        "selected_event.fraction." + event_type_name + ".by_equivalent_index";
    std::string desc =
        "Selected event count, for all events of a single event type by "
        "equivalent index. In the set of symmetrically equivalent events, "
        "events with the same equivalent index differ only by a translation. "
        "Requires selected_event." +
        event_type_name +
        ".by_equivalent_index selected event data is collected.";

    std::string histogram_name =
        "selected_event." + event_type_name + ".by_equivalent_index";
    auto const &event_f =
        _get_vector_int_histogram_function(calculation, name, histogram_name);
    bool sample_count = false;

    DiscreteVectorIntHistogramSamplingFunction hist_sampling_f(
        calculation, name, histogram_name, sample_count, *event_f.value_labels);

    f_list.emplace_back(name, desc, hist_sampling_f.component_names(),
                        hist_sampling_f.shape(), hist_sampling_f);
  }
  return f_list;
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM
