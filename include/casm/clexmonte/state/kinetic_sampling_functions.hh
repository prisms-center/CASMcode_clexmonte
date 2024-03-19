#ifndef CASM_clexmonte_state_kinetic_sampling_functions
#define CASM_clexmonte_state_kinetic_sampling_functions

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

/// \brief Make center of mass isotropic squared displacement sampling function
///     ("mean_R_squared_collective_isotropic")
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_collective_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make center of mass anisotropic squared displacement sampling
/// function
/// ("mean_R_squared_collective_anisotropic")
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_collective_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make tracer isotropic squared displacement sampling function
///     ("mean_R_squared_individual_isotropic")
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_individual_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make tracer anisotropic squared displacement sampling function
///     ("mean_R_squared_individual_anisotropic")
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_individual_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make isotropic Onsager kinetic coefficient sampling function
///     ("L_isotropic")
template <typename CalculationType>
state_sampling_function_type make_L_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make anisotropic Onsager kinetic coefficient sampling function
///     ("L_anisotropic")
template <typename CalculationType>
state_sampling_function_type make_L_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make isotropic tracer diffusion coefficient sampling function
///     ("D_tracer_isotropic")
template <typename CalculationType>
state_sampling_function_type make_D_tracer_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make anisotropic tracer diffusion coefficient sampling function
///     ("D_tracer_anisotropic")
template <typename CalculationType>
state_sampling_function_type make_D_tracer_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make delta_n_jumps(i) / n_atoms(i) ("jumps_per_atom_by_type")
template <typename CalculationType>
state_sampling_function_type make_jumps_per_atom_by_type_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make delta_n_jumps(i) / delta_n_events ("jumps_per_event_by_type")
template <typename CalculationType>
state_sampling_function_type make_jumps_per_event_by_type_f(
    std::shared_ptr<CalculationType> const &calculation);

/// \brief Make delta_n_jumps(i) / n_atoms(i) / delta_n_events
/// ("jumps_per_atom_per_event_by_type")
template <typename CalculationType>
state_sampling_function_type make_jumps_per_atom_per_event_by_type_f(
    std::shared_ptr<CalculationType> const &calculation);

// --- Inline definitions ---

/// \brief Make center of mass isotropic squared displacement sampling function
///     ("mean_R_squared_collective_isotropic")
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_collective_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_collective_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_individual_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_mean_R_squared_individual_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_L_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_L_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_D_tracer_isotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_D_tracer_anisotropic_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &kmc_data = calculation->kmc_data;
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
template <typename CalculationType>
state_sampling_function_type make_jumps_per_atom_by_type_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list =
            calculation->kmc_data.atom_name_index_list;
        auto const &n_jumps = calculation->occ_location->current_atom_n_jumps();

        auto const &sampling_fixture = *calculation->kmc_data.sampling_fixture;
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
template <typename CalculationType>
state_sampling_function_type make_jumps_per_event_by_type_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list =
            calculation->kmc_data.atom_name_index_list;
        auto const &n_jumps = calculation->occ_location->current_atom_n_jumps();

        auto const &sampling_fixture = *calculation->kmc_data.sampling_fixture;
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
template <typename CalculationType>
state_sampling_function_type make_jumps_per_atom_per_event_by_type_f(
    std::shared_ptr<CalculationType> const &calculation) {
  // Construct component_names && shape
  auto const &system = *calculation->system;
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
        auto const &system = *calculation->system;
        auto event_system = get_event_system(system);

        auto const &name_list = event_system->atom_name_list;
        auto const &name_index_list =
            calculation->kmc_data.atom_name_index_list;
        auto const &n_jumps = calculation->occ_location->current_atom_n_jumps();

        auto const &sampling_fixture = *calculation->kmc_data.sampling_fixture;
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

}  // namespace clexmonte
}  // namespace CASM

#endif
