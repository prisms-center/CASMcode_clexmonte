#include "casm/clexmonte/system/OccSystem.hh"

#include "casm/clexmonte/clex/Configuration.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"

namespace CASM {
namespace clexmonte {

/// \brief Constructor
OccSystem::OccSystem(
    std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
    composition::CompositionConverter const &_composition_converter,
    ClexData const &_formation_energy_clex_data)
    : shared_prim(_shared_prim),
      global_dof_info(clexulator::make_global_dof_info(*shared_prim)),
      local_dof_info(clexulator::make_local_dof_info(*shared_prim)),
      composition_converter(_composition_converter),
      composition_calculator(composition_converter.components(),
                             xtal::allowed_molecule_names(*shared_prim)),
      formation_energy_clex_data(_formation_energy_clex_data) {}

/// \brief Constructor
OccSystemSupercellData::OccSystemSupercellData(
    ClexData const &_formation_energy_clex_data,
    Eigen::Matrix3l const &_transformation_matrix_to_super)
    : supercell_neighbor_list(std::make_shared<clexulator::SuperNeighborList>(
          _transformation_matrix_to_super,
          *(_formation_energy_clex_data.prim_neighbor_list))),
      formation_energy_clex(supercell_neighbor_list,
                            _formation_energy_clex_data.clexulator,
                            _formation_energy_clex_data.eci) {}

// --- The following are used to construct a common interface between "System"
// data, in this case OccSystem, and templated CASM::clexmonte methods such as
// sampling function factory methods ---

namespace {

/// \brief Helper to get OccSystemSupercellData,
///     constructing as necessary
OccSystemSupercellData &get_supercell_data(
    OccSystem &data, Eigen::Matrix3l const &transformation_matrix_to_super) {
  auto it = data.supercell_data.find(transformation_matrix_to_super);
  if (it == data.supercell_data.end()) {
    it = data.supercell_data
             .emplace(std::piecewise_construct,
                      std::forward_as_tuple(transformation_matrix_to_super),
                      std::forward_as_tuple(data.formation_energy_clex_data,
                                            transformation_matrix_to_super))
             .first;
  }
  return it->second;
}

/// \brief Helper to get OccSystemSupercellData,
///     constructing as necessary
OccSystemSupercellData &get_supercell_data(
    OccSystem &data, monte::State<Configuration> const &state) {
  auto const &configuration = state.configuration;
  auto const &T = get_transformation_matrix_to_super(configuration);
  return get_supercell_data(data, T);
}

}  // namespace

/// \brief Helper to get std::shared_ptr<xtal::BasicStructure const>
std::shared_ptr<xtal::BasicStructure const> const &get_shared_prim(
    OccSystem &data) {
  return data.shared_prim;
}

/// \brief Helper to get composition::CompositionConverter
composition::CompositionConverter const &get_composition_converter(
    OccSystem &data) {
  return data.composition_converter;
}

/// \brief Helper to get composition::CompositionCalculator
composition::CompositionCalculator const &get_composition_calculator(
    OccSystem &data) {
  return data.composition_calculator;
}

/// \brief Helper to make the default configuration in a supercell
Configuration make_default_configuration(
    OccSystem const &data,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return Configuration(transformation_matrix_to_super,
                       clexulator::make_default_config_dof_values(
                           data.shared_prim->basis().size(),
                           transformation_matrix_to_super.determinant(),
                           data.global_dof_info, data.local_dof_info));
}

/// \brief Convert configuration from standard basis to prim basis
Configuration from_standard_values(
    OccSystem const &data,
    Configuration const &configuration_in_standard_basis) {
  Eigen::Matrix3l const &T =
      configuration_in_standard_basis.transformation_matrix_to_super;
  return Configuration(T, clexulator::from_standard_values(
                              configuration_in_standard_basis.dof_values,
                              data.shared_prim->basis().size(), T.determinant(),
                              data.global_dof_info, data.local_dof_info));
}

/// \brief Convert configuration from prim basis to standard basis
Configuration to_standard_values(
    OccSystem const &data, Configuration const &configuration_in_prim_basis) {
  Eigen::Matrix3l const &T =
      configuration_in_prim_basis.transformation_matrix_to_super;
  return Configuration(T, clexulator::to_standard_values(
                              configuration_in_prim_basis.dof_values,
                              data.shared_prim->basis().size(), T.determinant(),
                              data.global_dof_info, data.local_dof_info));
}

/// \brief Helper to get the ClexData for formation energy
ClexData &get_formation_energy_clex_data(OccSystem &data) {
  return data.formation_energy_clex_data;
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular configuration, constructing as necessary
///
/// \relates OccSystem
clexulator::ClusterExpansion &get_formation_energy_clex(
    OccSystem &data, Configuration const &configuration) {
  auto &clex =
      get_supercell_data(data, configuration.transformation_matrix_to_super)
          .formation_energy_clex;
  clex.set(&configuration.dof_values);
  return clex;
}

/// \brief Helper to get the correct clexulator::ClusterExpansion for a
///     particular state, constructing as necessary
///
/// \relates OccSystem
clexulator::ClusterExpansion &get_formation_energy_clex(
    OccSystem &data, monte::State<Configuration> const &state) {
  auto &clex = get_supercell_data(data, state).formation_energy_clex;
  clex.set(&state.configuration.dof_values);
  return clex;
}

}  // namespace clexmonte
}  // namespace CASM
