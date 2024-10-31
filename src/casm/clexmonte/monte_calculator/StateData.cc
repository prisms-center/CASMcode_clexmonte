
#include "casm/clexmonte/monte_calculator/StateData.hh"

#include "casm/clexmonte/state/Configuration.hh"

namespace CASM {
namespace clexmonte {

StateData::StateData(std::shared_ptr<system_type> _system,
                     state_type const *_state,
                     monte::OccLocation const *_occ_location)
    : system(_system), state(_state), occ_location(_occ_location) {
  if (system == nullptr) {
    throw std::runtime_error("Error constructing StateData: system==nullptr");
  }
  if (state == nullptr) {
    throw std::runtime_error("Error constructing StateData: state==nullptr");
  }
  //  if (occ_location == nullptr) {
  //    throw std::runtime_error(
  //        "Error constructing StateData: occ_location==nullptr");
  //  }

  transformation_matrix_to_super = get_transformation_matrix_to_super(*state);
  n_unitcells = transformation_matrix_to_super.determinant();
  convert = &get_index_conversions(*system, *state);

  // make supercell_neighbor_list
  auto supercell_neighbor_list = get_supercell_neighbor_list(*system, *state);
  if (supercell_neighbor_list == nullptr) {
    throw std::runtime_error(
        "Error constructing StateData: empty supercell neighbor list");
  }

  // make corr
  for (auto const &pair : system->basis_sets) {
    auto const &key = pair.first;
    // create a copy of pair.second,
    // to ensure clexulator is pointing at this->state
    std::shared_ptr<clexulator::Clexulator> clexulator =
        std::make_shared<clexulator::Clexulator>(*pair.second);
    auto _corr = std::make_shared<clexulator::Correlations>(
        supercell_neighbor_list, clexulator);
    _corr->set(&get_dof_values(*state));
    corr.emplace(key, _corr);
  }

  // make local_corr
  for (auto const &pair : system->local_basis_sets) {
    auto const &key = pair.first;
    // create a copy of pair.second,
    // to ensure clexulator are pointing at this->state
    std::shared_ptr<std::vector<clexulator::Clexulator>> _local_clexulator =
        std::make_shared<std::vector<clexulator::Clexulator>>();
    for (clexulator::Clexulator const &_equiv : *pair.second) {
      _local_clexulator->push_back(_equiv);
    }
    auto _local_corr = std::make_shared<clexulator::LocalCorrelations>(
        supercell_neighbor_list, _local_clexulator);
    _local_corr->set(&get_dof_values(*state));
    local_corr.emplace(key, _local_corr);
  }

  // make clex
  for (auto const &pair : system->clex_data) {
    auto const &key = pair.first;
    auto const &data = pair.second;

    // create a copy of *get_basis_set(system, data.basis_set_name),
    // to ensure clexulator is pointing at this->state
    auto _tmp = get_basis_set(*system, data.basis_set_name);
    auto _clexulator = std::make_shared<clexulator::Clexulator>(*_tmp);

    // construct ClusterExpansion
    auto _clex = std::make_shared<clexulator::ClusterExpansion>(
        supercell_neighbor_list, _clexulator, data.coefficients);
    set(*_clex, *state);
    clex.emplace(key, _clex);
  }

  // make multiclex
  for (auto const &pair : system->multiclex_data) {
    auto const &key = pair.first;
    auto const &data = pair.second;

    // create a copy of *get_basis_set(system, data.basis_set_name),
    // to ensure clexulator is pointing at this->state
    auto _tmp = get_basis_set(*system, data.basis_set_name);
    auto _clexulator = std::make_shared<clexulator::Clexulator>(*_tmp);

    // construct MultiClusterExpansion
    auto _multiclex = std::make_shared<clexulator::MultiClusterExpansion>(
        supercell_neighbor_list, _clexulator, data.coefficients);
    set(*_multiclex, *state);
    auto _glossary = data.coefficients_glossary;
    multiclex.emplace(key, std::make_pair(_multiclex, _glossary));
  }

  // make local_clex
  for (auto const &pair : system->local_clex_data) {
    auto const &key = pair.first;
    auto const &data = pair.second;

    // create a copy of *get_local_basis_set(system, data.basis_set_name),
    // to ensure clexulator are pointing at this->state
    auto _tmp = get_local_basis_set(*system, data.local_basis_set_name);
    std::shared_ptr<std::vector<clexulator::Clexulator>> _local_clexulator =
        std::make_shared<std::vector<clexulator::Clexulator>>();
    for (clexulator::Clexulator const &_equiv : *_tmp) {
      _local_clexulator->push_back(_equiv);
    }

    // construct LocalClusterExpansion
    auto _local_clex = std::make_shared<clexulator::LocalClusterExpansion>(
        supercell_neighbor_list, _local_clexulator, data.coefficients);
    set(*_local_clex, *state);
    local_clex.emplace(key, _local_clex);
  }

  // make local_multiclex
  for (auto const &pair : system->local_multiclex_data) {
    auto const &key = pair.first;
    auto const &data = pair.second;
    // create a copy of *get_local_basis_set(system, data.basis_set_name),
    // to ensure clexulator are pointing at this->state
    auto _tmp = get_local_basis_set(*system, data.local_basis_set_name);
    std::shared_ptr<std::vector<clexulator::Clexulator>> _local_clexulator =
        std::make_shared<std::vector<clexulator::Clexulator>>();
    for (clexulator::Clexulator const &_equiv : *_tmp) {
      _local_clexulator->push_back(_equiv);
    }

    // construct MultiLocalClusterExpansion
    auto _local_multiclex =
        std::make_shared<clexulator::MultiLocalClusterExpansion>(
            supercell_neighbor_list, _local_clexulator, data.coefficients);
    set(*_local_multiclex, *state);
    auto _glossary = data.coefficients_glossary;
    local_multiclex.emplace(key, std::make_pair(_local_multiclex, _glossary));
  }

  // make order_parameters
  for (auto const &pair : system->dof_spaces) {
    auto const &key = pair.first;
    auto const &definition = *pair.second;
    auto _order_parameter =
        std::make_shared<clexulator::OrderParameter>(definition);
    _order_parameter->update(convert->transformation_matrix_to_super(),
                             convert->index_converter(),
                             &get_dof_values(*state));
    order_parameters.emplace(key, _order_parameter);
  }

  // make local orbit composition calculators
  for (auto const &pair : system->local_orbit_composition_calculator_data) {
    auto const &key = pair.first;
    auto const &data = pair.second;

    auto const &composition_calculator = get_composition_calculator(*system);
    auto const &orbits =
        get_local_basis_set_cluster_info(*system, data->local_basis_set_name)
            ->orbits;
    auto prim_nlist = system->prim_neighbor_list;
    auto supercell_nlist = get_supercell_neighbor_list(*system, *state);
    auto const &supercell_index_converter =
        get_index_conversions(*system, *state).index_converter();
    clexulator::ConfigDoFValues const *dof_values = &get_dof_values(*state);

    auto _local_orbit_composition_calculator =
        std::make_shared<LocalOrbitCompositionCalculator>(
            orbits, data->orbits_to_calculate, data->combine_orbits, prim_nlist,
            supercell_nlist, supercell_index_converter, composition_calculator,
            dof_values);
    local_orbit_composition_calculators.emplace(
        key, _local_orbit_composition_calculator);
  }
}

}  // namespace clexmonte
}  // namespace CASM
