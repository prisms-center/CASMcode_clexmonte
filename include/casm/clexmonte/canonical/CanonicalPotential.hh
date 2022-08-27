#ifndef CASM_clexmonte_canonical_CanonicalPotential
#define CASM_clexmonte_canonical_CanonicalPotential

#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexulator/ClusterExpansion.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Implements potential for canonical Monte Carlo
class CanonicalPotential {
 public:
  CanonicalPotential(
      std::shared_ptr<clexulator::ClusterExpansion> _formation_energy_clex)
      : m_formation_energy_clex(_formation_energy_clex) {}

  /// \brief Reset pointer to state currently being calculated
  void set(monte::State<Configuration> const *state) {
    m_state = state;
    clexmonte::set(*m_formation_energy_clex, *m_state);
  }

  /// \brief Pointer to state currently being calculated
  monte::State<Configuration> const *get() const { return m_state; }

  /// \brief Calculate (extensive) cluster expansion value
  double extensive_value() {
    return m_formation_energy_clex->extensive_value();
  }

  /// \brief Calculate change in (extensive) cluster expansion value due to a
  ///     series of occupation changes
  double occ_delta_extensive_value(std::vector<Index> const &linear_site_index,
                                   std::vector<int> const &new_occ) {
    return m_formation_energy_clex->occ_delta_value(linear_site_index, new_occ);
  }

 private:
  /// State to use
  monte::State<Configuration> const *m_state;

  /// Formation energy cluster expansion calculator;
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
};

/// \brief Set potential calculator so it evaluates using `state`
inline void set(CanonicalPotential &potential,
                monte::State<Configuration> const &state) {
  potential.set(&state);
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
