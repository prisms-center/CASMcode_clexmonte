#ifndef CASM_clexmonte_kmc_PrimEventCalculator
#define CASM_clexmonte_kmc_PrimEventCalculator

#include <memory>

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/kmc/conditions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace kmc {

/// \brief Event rate calculation for a particular KMC event
class PrimEventCalculator {
 public:
  PrimEventCalculator(
      std::shared_ptr<clexulator::ClusterExpansion> formation_energy_clex,
      std::shared_ptr<clexulator::MultiLocalClusterExpansion> event_clex)
      : m_formation_energy_clex(formation_energy_clex),
        m_event_clex(event_clex) {}

  /// \brief Calculate the state of an event
  ///
  /// \param state Stores whether the event is allowed, is "normal",
  ///     energy barriers, and event rate
  /// \param conditions Holds any conditions (i.e. beta) necessary for
  ///     calculating event state
  /// \param event_data Holds information about the particular translational
  ///     instance of the event, such as linear sites indices and linear
  ///     unitcell index, necessary for calculating event state.
  /// \param event_data Holds information about the event that does not
  ///     depend on the particular translational instance, such as the
  ///     initial and final occupation variables.
  void calculate_event_state(EventState &state, Conditions const &conditions,
                             EventData const &event_data,
                             PrimEventData const &prim_event_data) const {
    clexulator::ConfigDoFValues const *dof_values =
        m_formation_energy_clex->get();

    int i = 0;
    for (Index l : event_data.event.linear_site_index) {
      if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
        state.is_allowed = false;
        state.rate = 0.0;
        return;
      }
      ++i;
    }
    state.is_allowed = true;

    // calculate change in energy to final state
    state.dE_final = m_formation_energy_clex->occ_delta_value(
        event_data.event.linear_site_index, prim_event_data.occ_final);

    // calculate KRA and attempt frequency
    Eigen::VectorXd const &event_values = m_event_clex->values(
        event_data.unitcell_index, prim_event_data.equivalent_index);
    state.Ekra =
        event_values[static_cast<unsigned long>(EVENT_CLEX_INDEX::KRA)];
    state.freq =
        event_values[static_cast<unsigned long>(EVENT_CLEX_INDEX::FREQ)];

    // calculate energy in activated state, check if "normal", calculate rate
    state.dE_activated = state.dE_final * 0.5 + state.Ekra;
    state.is_normal =
        (state.dE_activated > 0.0) && (state.dE_activated > state.dE_final);
    if (state.dE_activated < state.dE_final)
      state.dE_activated = state.dE_final;
    if (state.dE_activated < 0.0) state.dE_activated = 0.0;
    state.rate = state.freq * exp(-conditions.beta * state.dE_activated);
  }

 private:
  // --- reference to current supercell ---
  std::shared_ptr<clexulator::ClusterExpansion> m_formation_energy_clex;
  std::shared_ptr<clexulator::MultiLocalClusterExpansion> m_event_clex;
};

/// \brief Construct PrimEventCalculator list
template <typename SystemType>
std::vector<kmc::PrimEventCalculator> make_prim_event_calculators(
    SystemType &system, monte::State<clexmonte::Configuration> const &state,
    std::vector<PrimEventData> const &prim_event_list) {
  std::vector<kmc::PrimEventCalculator> prim_event_calculators;
  for (auto const &prim_event_data : prim_event_list) {
    prim_event_calculators.emplace_back(
        get_clex(system, state, "formation_energy"),
        get_local_multiclex(system, state, prim_event_data.event_type_name));
  }
  return prim_event_calculators;
}

}  // namespace kmc
}  // namespace clexmonte
}  // namespace CASM

#endif
