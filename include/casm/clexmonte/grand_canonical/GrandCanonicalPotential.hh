#ifndef CASM_clexmonte_grandcanonical_GrandCanonicalPotential
#define CASM_clexmonte_grandcanonical_GrandCanonicalPotential

#include "casm/clexmonte/Configuration.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/monte/Conversions.hh"

namespace CASM {
namespace clexmonte {
namespace grandcanonical {

/// \brief Implements potential for (semi) grand canonical Monte Carlo
///
/// Notes:
/// - Expects Monte Carlo states to have "chem_pot" or "param_chem_pot"
///   conditions
class GrandCanonicalPotential {
  GrandCanonicalPotential(
      std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex,
      composition::CompositionConverter const &composition_converter,
      composition::CompositionCalculator const &composition_calculator,
      monte::Conversions const &convert)
      : m_formation_energy_clex(formation_energy_clex),
        m_composition_calculator(composition_calculator),
        m_composition_calculator(composition_calculator),
        m_convert(convert) {}

  /// \brief Reset pointer to state currently being calculated
  void set(monte::State<Configuration> const *state) {
    m_state = state;
    m_n_unitcells = get_transformation_matrix_to_super(m_state->configuration)
                        .determinant();

    // get parametric chemical potential
    if (state.conditions.scalar_values.count("param_chem_pot")) {
      m_param_chem_pot = state.conditions.scalar_values.at("param_chem_pot");
    } else {
      throw std::runtime_error(
          "Error in clexmonte::GrandCanonicalPotential: missing "
          "\"param_chem_pot\" conditions.");
    }
    if (m_param_chem_pot.size() !=
        m_composition_converter.independent_compositions()) {
      throw std::runtime_error(
          "Error in clexmonte::GrandCanonicalPotential: mismatch between "
          "\"param_chem_pot\" size and composition axes");
    }

    // Make matrix of Mij = mu_i - mu_j
    m_exchange_chem_pot = make_exchange_chemical_potential(
        m_param_chem_pot, m_composition_converter);
  }

  /// \brief Pointer to state currently being calculated
  monte::State<Configuration> const *get() const { return m_configuration; }

  /// \brief Calculate (extensive) cluster expansion value
  double extensive_value() {
    auto const &occupation = get_occupation(m_state->configuration);
    Eigen::VectorXd comp_n =
        m_composition_calculator.mean_num_each_component(occupation);
    Eigen::VectorXd comp_x = m_composition_converter.param_composition(comp_n);
    return m_formation_energy_clex.extensive_value() -
           m_n_unitcells * m_param_chem_pot.dot(comp_x);
  }

  /// \brief Calculate change in (extensive) cluster expansion value due to a
  ///     single occupation change
  double occ_delta_extensive_value(Index linear_site_index, int new_occ) {
    Index asym_index = m_convert.l_to_asym(linear_site_index);
    Index new_species = m_convert.species_index(asym_index, new_occ);
    Index curr_occ = get_occupation(m_state->configuration)(linear_site_index);
    Index curr_species = m_convert.species_index(asym_index, curr_occ);
    return m_formation_energy_clex.occ_delta_value(linear_site_index, new_occ) -
           m_exchange_chem_pot(new_species, curr_species);
  }

  /// \brief Calculate change in (extensive) cluster expansion value due to a
  ///     series of occupation changes
  double occ_delta_extensive_value(std::vector<Index> const &linear_site_index,
                                   std::vector<int> const &new_occ) {
    double dE =
        m_formation_energy_clex.occ_delta_value(linear_site_index, new_occ);
    Eigen::VectorXi const &occupation = get_occupation(m_state->configuration);
    for (Index i = 0; i < linear_site_index.size(); ++i) {
      Index l = linear_site_index[i];
      Index asym_index = m_convert.l_to_asym(l);
      Index new_species = m_convert.species_index(asym_index, new_occ[i]);
      Index curr_occ = occupation(l);
      Index curr_species = m_convert.species_index(asym_index, curr_occ);
      dE -= m_exchange_chem_pot(new_species, curr_species);
    }
    return dE;
  }

 private:
  /// Current state to use in calculations
  monte::State<ConfigType> const *m_state;

  /// Number of unit cells, from current state
  double m_n_unitcells;

  /// Parametric chemical potential, from current state
  Eigen::VectorXd m_param_chem_pot;

  /// Matrix, \f$M\f$, with values \f$M(i,j) = \mu_i - \mu_j\f$, calculated
  /// from `m_param_chem_pot`.
  Eigen::MatrixXd m_exchange_chem_pot;

  /// Formation energy cluster expansion calculator;
  clexulator::ClusterExpansion m_formation_energy_clex;

  /// Number per unit cell / parametric composition conversions
  composition::CompositionConverter m_composition_converter;

  /// Composition calculator
  composition::CompositionCalculator m_composition_calculator;

  /// Performs index conversions
  monte::Conversions m_convert;
};

}  // namespace grandcanonical
}  // namespace clexmonte
}  // namespace CASM

#endif
