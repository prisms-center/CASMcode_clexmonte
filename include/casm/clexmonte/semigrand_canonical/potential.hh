#ifndef CASM_clexmonte_semigrand_canonical_potential
#define CASM_clexmonte_semigrand_canonical_potential

#include <random>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexulator/ClusterExpansion.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

/// \brief Implements potential for semi-grand canonical Monte Carlo
class SemiGrandCanonicalPotential {
 public:
  SemiGrandCanonicalPotential(std::shared_ptr<system_type> _system);

  /// \brief Pointer to configuration currently being calculated
  clexulator::ConfigDoFValues const *get() const;

  /// \brief Reset pointer to state currently being calculated
  void set(state_type const *state,
           std::shared_ptr<SemiGrandCanonicalConditions> conditions);

  /// \brief Pointer to current state
  state_type const *state() const;

  /// \brief Pointer to current conditions, valid for current state
  std::shared_ptr<SemiGrandCanonicalConditions> const &conditions() const;

  /// \brief Pointer to formation energy cluster expansion calculator, valid for
  /// current state
  std::shared_ptr<clexulator::ClusterExpansion> const &formation_energy() const;

  /// \brief Calculate (per_supercell) semi-grand potential value
  double per_supercell();

  /// \brief Calculate (per_unitcell) semi-grand potential value
  double per_unitcell();

  /// \brief Calculate change in (per_supercell) semi-grand potential value due
  ///     to a series of occupation changes
  double occ_delta_per_supercell(std::vector<Index> const &linear_site_index,
                                 std::vector<int> const &new_occ);

 private:
  /// System pointer
  std::shared_ptr<system_type> m_system;

  /// State to use
  state_type const *m_state;

  /// Formation energy cluster expansion calculator, depends on current state
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;

  /// Number of unit cells, depends on current state
  double m_n_unitcells;

  /// Conditions, depends on current state
  std::shared_ptr<SemiGrandCanonicalConditions> m_conditions;

  /// Index conversions, depends on current state
  monte::Conversions const *m_convert;
};

typedef SemiGrandCanonicalPotential potential_type;

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
