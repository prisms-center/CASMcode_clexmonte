
#include "casm/clexmonte/monte_calculator/sampling_functions.hh"

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/misc/to_json.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/clexulator/Correlations.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/sampling/StateSamplingFunction.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {
namespace monte_calculator {

/// \brief Make temperature sampling function ("temperature")
///
/// Requires:
/// - "temperature" is a scalar state condition
state_sampling_function_type make_temperature_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return state_sampling_function_type(
      "temperature", "Temperature (K)", {},  // scalar,
      [calculation]() {
        return scalar_conditions(calculation, "temperature");
      });
}

/// \brief Make mol composition sampling function ("mol_composition")
state_sampling_function_type make_mol_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  auto const &system = get_system(calculation);
  auto const &components = get_composition_converter(system).components();
  std::vector<Index> shape;
  shape.push_back(components.size());
  return state_sampling_function_type(
      "mol_composition",
      "Number of each component (normalized per primitive cell)",
      components,  // component names
      shape, [calculation]() {
        auto const &system = get_system(calculation);
        auto const &state = get_state(calculation);
        Eigen::VectorXi const &occupation = get_occupation(state);
        return get_composition_calculator(system).mean_num_each_component(
            occupation);
      });
}

/// \brief Make parametric composition sampling function ("param_composition")
state_sampling_function_type make_param_composition_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  auto const &system = get_system(calculation);
  // name param_composition components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(system);
  std::vector<std::string> component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    component_names.push_back(composition_converter.comp_var(i));
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "param_composition", "Parametric composition",
      component_names,  // component names
      shape, [calculation]() {
        auto const &system = get_system(calculation);
        auto const &state = get_state(calculation);
        composition::CompositionCalculator const &composition_calculator =
            get_composition_calculator(system);
        composition::CompositionConverter const &composition_converter =
            get_composition_converter(system);

        Eigen::VectorXi const &occupation = get_occupation(state);
        Eigen::VectorXd mol_composition =
            composition_calculator.mean_num_each_component(occupation);
        return composition_converter.param_composition(mol_composition);
      });
}

/// \brief Make parametric chemical potential sampling function
/// ("param_chem_pot")
///
/// Requires:
/// - "param_chem_pot" is a vector state condition
state_sampling_function_type make_param_chem_pot_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  auto const &system = get_system(calculation);
  // name param_chem_pot components "a", "b", ... for each independent
  // composition axis
  composition::CompositionConverter const &composition_converter =
      get_composition_converter(system);
  std::vector<std::string> component_names;
  for (Index i = 0; i < composition_converter.independent_compositions(); ++i) {
    component_names.push_back(composition_converter.comp_var(i));
  }
  std::vector<Index> shape;
  shape.push_back(component_names.size());

  return state_sampling_function_type(
      "param_chem_pot",
      "Chemical potential conjugate to parametric composition axes",
      component_names,  // component names
      shape, [calculation]() {
        return vector_conditions(calculation, "param_chem_pot");
      });
}

/// \brief Make correlations sampling function ("corr.<key>")
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::corr, a basis set name
state_sampling_function_type make_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  std::vector<Index> shape;
  Index size = get_basis_set(get_system(calculation), key)->corr_size();
  shape.push_back(size);

  return state_sampling_function_type(
      std::string("corr.") + key,
      "Correlations values (normalized per primitive cell)", shape,
      [calculation, key]() {
        auto &correlations = calculation->state_data()->corr.at(key);
        auto const &per_supercell_corr = correlations->per_supercell();
        return correlations->per_unitcell(per_supercell_corr);
      });
}

/// \brief Make cluster expansion value sampling function ("clex.<key>")
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::clex, a cluster expansion name
state_sampling_function_type make_clex_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  return state_sampling_function_type(
      std::string("clex.") + key,
      "Cluster expansion value (normalized per primitive cell)", {},  // scalar
      [calculation, key]() {
        Eigen::VectorXd value(1);
        value(0) = calculation->state_data()->clex.at(key)->per_unitcell();
        return value;
      });
}

/// \brief Make multi-cluster expansion value sampling function
/// ("multiclex.<key>")
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::multiclex, a multi-cluster expansion name
state_sampling_function_type make_multiclex_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  auto const &data = get_multiclex_data(get_system(calculation), key);
  std::vector<Index> shape;
  Index size = data.coefficients.size();
  shape.push_back(size);

  std::vector<std::string> component_names;
  for (Index i = 0; i < size; ++i) {
    component_names.push_back(std::to_string(i));
  }
  for (auto const &pair : data.coefficients_glossary) {
    component_names[pair.second] = pair.first;
  }

  return state_sampling_function_type(
      std::string("clex.") + key,
      "Mulit-cluster expansion value (normalized per primitive cell)",
      component_names,
      shape,  // scalar
      [calculation, key]() {
        return calculation->state_data()
            ->multiclex.at(key)
            .first->per_unitcell();
      });
}

/// \brief Make non-zero coefficients correlations sampling function
///     ("clex.<key>.sparse_corr")
state_sampling_function_type make_clex_sparse_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  auto const &system = get_system(calculation);
  // get sparse correlations size
  auto const &clex_data = get_clex_data(system, key);
  auto const &coeff = clex_data.coefficients;
  std::vector<Index> shape;
  shape.push_back(coeff.index.size());
  std::vector<std::string> component_names;
  for (unsigned int index : coeff.index) {
    component_names.push_back(std::to_string(index));
  }

  return state_sampling_function_type(
      std::string("clex.") + key + std::string(".sparse_corr"),
      "Cluster expansion correlations, for non-zero coefficients (normalized "
      "per primitive cell)",
      component_names, shape, [calculation, key]() {
        auto &correlations =
            calculation->state_data()->clex.at(key)->correlations();
        auto const &per_supercell_corr = correlations.per_supercell();
        Eigen::VectorXd all_corr =
            correlations.per_unitcell(per_supercell_corr);
        auto const &indices = correlations.correlation_indices();
        Eigen::VectorXd sparse_corr(indices.size());
        Index i = 0;
        for (unsigned int index : indices) {
          sparse_corr(i) = all_corr(index);
          ++i;
        }
        return sparse_corr;
      });
}

/// \brief Make non-zero coefficients correlations sampling function
///     ("multiclex.<key>.sparse_corr")
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::multiclex, a multi-cluster expansion name
state_sampling_function_type make_multiclex_sparse_corr_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  auto const &system = get_system(calculation);
  // get combined sparse correlations size
  auto const &multiclex_data = get_multiclex_data(system, key);
  std::set<unsigned int> combined_index;
  for (auto const &coeff : multiclex_data.coefficients) {
    for (auto index : coeff.index) {
      combined_index.insert(index);
    }
  }
  std::vector<Index> shape;
  shape.push_back(combined_index.size());
  std::vector<std::string> component_names;
  for (unsigned int index : combined_index) {
    component_names.push_back(std::to_string(index));
  }

  return state_sampling_function_type(
      std::string("multiclex.") + key + std::string(".sparse_corr"),
      "Cluster expansion correlations, for non-zero coefficients (normalized "
      "per primitive cell)",
      component_names, shape, [calculation, key]() {
        auto &correlations =
            calculation->state_data()->multiclex.at(key).first->correlations();
        auto const &per_supercell_corr = correlations.per_supercell();
        Eigen::VectorXd all_corr =
            correlations.per_unitcell(per_supercell_corr);
        auto const &indices = correlations.correlation_indices();
        Eigen::VectorXd sparse_corr(indices.size());
        Index i = 0;
        for (unsigned int index : indices) {
          sparse_corr(i) = all_corr(index);
          ++i;
        }
        return sparse_corr;
      });
}

/// \brief Make potential energy sampling function (allows user-specified
/// <label>)
///
/// Notes:
/// - This uses calculation->potential->per_unitcell()
///
/// \param calculation Monte Carlo calculator
/// \param label Name to give the sampling function
/// \param desc Description to give the sampling function
state_sampling_function_type make_potential_energy_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string label,
    std::string desc) {
  // desc = "Potential energy of the state (normalized per primitive cell)",
  return state_sampling_function_type(
      label, desc, {},  // scalar
      [calculation]() {
        Eigen::VectorXd value(1);
        value(0) = calculation->potential().per_unitcell();
        return value;
      });
}

/// \brief Make order parameter sampling function ("order_parameter.<key>")
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::dof_spaces
state_sampling_function_type make_order_parameter_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  clexulator::DoFSpace const &dof_space =
      *get_system(calculation).dof_spaces.at(key);
  std::string name = "order_parameter." + key;
  std::string desc = "Order parameters";

  return state_sampling_function_type(
      name, desc, {dof_space.subspace_dim},  // vector size
      [calculation, key]() {
        return calculation->state_data()->order_parameters.at(key)->value();
      });
}

/// \brief Make order parameter magnitudes by subspace sampling function
/// ("order_parameter.<key>.subspace_magnitudes")
///
/// Creates a "order_parameter.<key>.subspace_magnitudes" function for the
/// specified DoFSpace, using the subspaces specified in the
/// `calculation->system()->dof_subspaces` map.
///
/// Example system input JSON, to measure the magnitude of the order parameter
/// in four distinct subspaces:
/// \code
/// {
///   ...
///   "dof_spaces": {
///     "0": "dof_spaces/dof_space.0.json"
///   },
///   "dof_subspaces": {
///     "0": [
///       [0],  # subspace formed by DoFSpace basis vector 0
///       [1, 2], # subspace formed by DoFSpace basis vectors 1, 2
///       [3, 4, 5], # subspace formed by DoFSpace basis vector 3, 4, 5
///       [6, 7, 8, 9, 10, 11]  # subspace formed by DoFSpace basis vector 6-11
///     ]
///   }
///   ...
/// }
/// \endcode
///
/// \param calculation Monte Carlo calculator
/// \param key Key into StateData::dof_spaces and StateData::dof_subspaces
state_sampling_function_type make_subspace_order_parameter_f(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  std::string name = "order_parameter." + key + ".subspace_magnitudes";
  std::string desc = "Order parameter magnitudes by subspace";

  auto const &system = get_system(calculation);
  auto it = system.dof_subspaces.find(key);
  if (it == system.dof_subspaces.end()) {
    std::stringstream msg;
    msg << "Error: no dof_subspaces for '" << key << "'";
    throw std::runtime_error(msg.str());
  }
  Index n_subspaces = it->second.size();
  return state_sampling_function_type(
      name, desc, {n_subspaces},  // vector size
      [calculation, key]() {
        Eigen::VectorXd eta =
            calculation->state_data()->order_parameters.at(key)->value();

        auto const &system = get_system(calculation);
        auto const &subspaces = system.dof_subspaces.at(key);
        Eigen::VectorXd eta_subspace = Eigen::VectorXd::Zero(subspaces.size());
        for (int i = 0; i < subspaces.size(); ++i) {
          double x = 0.0;
          for (int j : subspaces[i]) {
            if (j < 0 || j >= eta.size()) {
              throw std::runtime_error("Invalid dof_subspaces");
            }
            x += eta(j) * eta(j);
          }
          eta_subspace(i) = sqrt(x);
        }
        return eta_subspace;
      });
}

/// \brief Make configuration sampling function
json_state_sampling_function_type make_config_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::string name = "config";
  std::string desc = "The Monte Carlo configuration as JSON";
  return json_state_sampling_function_type(
      name, desc, [calculation]() -> jsonParser {
        auto const &state = get_state(calculation);
        return qto_json(state.configuration);
      });
}

/// \brief Make sampling functions that all methods can use
std::vector<state_sampling_function_type> common_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation,
    std::string potential_energy_name, std::string potential_energy_desc) {
  std::vector<state_sampling_function_type> functions = {
      make_temperature_f(calculation), make_mol_composition_f(calculation),
      make_param_composition_f(calculation)};

  auto const &system = *calculation->system();
  // make_corr_f,
  for (auto const &pair : system.basis_sets) {
    functions.push_back(make_corr_f(calculation, pair.first));
  }
  // make_clex_f, make_clex_sparse_corr_f,
  for (auto const &pair : system.clex_data) {
    functions.push_back(make_clex_f(calculation, pair.first));
    functions.push_back(make_clex_sparse_corr_f(calculation, pair.first));
  }
  // make_multiclex_f, make_multiclex_sparse_corr_f
  for (auto const &pair : system.multiclex_data) {
    functions.push_back(make_multiclex_f(calculation, pair.first));
    functions.push_back(make_multiclex_sparse_corr_f(calculation, pair.first));
  }
  // make_order_parameter_f
  for (auto const &pair : system.dof_spaces) {
    functions.push_back(make_order_parameter_f(calculation, pair.first));
  }
  // make_subspace_order_parameter_f
  for (auto const &pair : system.dof_subspaces) {
    if (!system.dof_spaces.count(pair.first)) {
      std::stringstream msg;
      msg << "Error: dof_subspaces includes '" << pair.first
          << "' which is not found in dof_spaces.";
      throw std::runtime_error(msg.str());
    }
    functions.push_back(
        make_subspace_order_parameter_f(calculation, pair.first));
  }
  // potential_energy / generalized_enthalpy / canonical_energy /
  // semigrand_canonical_energy / etc.
  functions.push_back(monte_calculator::make_potential_energy_f(
      calculation, potential_energy_name, potential_energy_desc));
  return functions;
}

/// \brief Make json sampling functions that all methods can use
std::vector<json_state_sampling_function_type> common_json_sampling_functions(
    std::shared_ptr<MonteCalculator> const &calculation) {
  std::vector<json_state_sampling_function_type> functions = {
      make_config_f(calculation)};
  return functions;
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM
