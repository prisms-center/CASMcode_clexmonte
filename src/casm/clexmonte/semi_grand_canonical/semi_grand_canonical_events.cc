
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical_events.hh"

#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

CompleteEventCalculator::CompleteEventCalculator(
    std::shared_ptr<SemiGrandCanonicalPotential> _potential,
    std::vector<PrimEventData> const &_prim_event_list,
    std::map<EventID, EventData> const &_event_list)
    : prim_event_list(_prim_event_list),
      event_list(_event_list),
      potential(_potential) {}

/// \brief Get CASM::monte::OccEvent corresponding to given event ID
double CompleteEventCalculator::calculate_rate(EventID const &id) {
  EventData const &event_data = event_list.at(id);
  PrimEventData const &prim_event_data =
      prim_event_list.at(id.prim_event_index);

  /// ---

  clexulator::ConfigDoFValues const *dof_values = potential->get();

  int i = 0;
  for (Index l : event_data.event.linear_site_index) {
    if (dof_values->occupation(l) != prim_event_data.occ_init[i]) {
      event_state.is_allowed = false;
      event_state.rate = 0.0;
      return event_state.rate;
    }
    ++i;
  }
  event_state.is_allowed = true;

  // calculate change in energy to final state
  event_state.dE_final = potential->occ_delta_extensive_value(
      event_data.event.linear_site_index, prim_event_data.occ_final);

  // calculate rate
  if (event_state.dE_final <= 0.0) {
    event_state.rate = 1.0;
  } else {
    event_state.rate =
        exp(-potential->conditions()->beta * event_state.dE_final);
  }

  /// ---

  return event_state.rate;
}

namespace {

/// \brief Convert monte::OccCanditate to occ_events::OccPosition
occ_events::OccPosition _make_occ_position(
    occ_events::OccSystem const &event_system,
    monte::Conversions const &convert, monte::OccCandidate const &cand) {
  Index b = *convert.asym_to_b(cand.asym).begin();
  return event_system.make_molecule_position(
      xtal::UnitCellCoord(b, 0, 0, 0),
      convert.occ_index(cand.asym, cand.species_index));
}

/// \brief Convert monte::OccSwap to occ_events::OccEvent
occ_events::OccEvent _make_grand_canonical_swap_event(
    monte::OccSwap const &swap, occ_events::OccSystem const &event_system,
    monte::Conversions const &convert) {
  occ_events::OccPosition pos_a =
      _make_occ_position(event_system, convert, swap.cand_a);
  occ_events::OccPosition resevoir_a =
      event_system.make_molecule_in_resevoir_position(
          event_system.get_chemical_index(pos_a));
  occ_events::OccPosition pos_b =
      _make_occ_position(event_system, convert, swap.cand_b);
  occ_events::OccPosition resevoir_b =
      event_system.make_molecule_in_resevoir_position(
          event_system.get_chemical_index(pos_b));

  occ_events::OccTrajectory traj_a({pos_a, resevoir_a});
  occ_events::OccTrajectory traj_b({resevoir_b, pos_b});

  return occ_events::OccEvent({traj_a, traj_b});
}

/// \brief Make OccEventTypeData for semi-grand canonical events
std::map<std::string, OccEventTypeData> _make_event_type_data(
    std::shared_ptr<system_type> system, state_type const &state,
    std::vector<monte::OccSwap> const &grand_canonical_swaps) {
  auto event_system = get_event_system(*system);
  monte::Conversions const &convert = get_index_conversions(*system, state);

  std::map<std::string, OccEventTypeData> event_type_data;
  auto const &occevent_symgroup_rep = get_occevent_symgroup_rep(*system);
  for (monte::OccSwap const &swap : grand_canonical_swaps) {
    // do not repeat forward and reverse -
    //   reverse will be constructed by make_prim_event_list(
    if (swap.cand_a.species_index < swap.cand_b.species_index) {
      continue;
    }
    occ_events::OccEvent event =
        _make_grand_canonical_swap_event(swap, *event_system, convert);
    std::set<occ_events::OccEvent> orbit =
        occ_events::make_prim_periodic_orbit(event, occevent_symgroup_rep);

    std::string event_type_name =
        "swap-" + std::to_string(swap.cand_a.species_index) + "-" +
        std::to_string(swap.cand_b.species_index);
    event_type_data[event_type_name].events =
        std::vector<occ_events::OccEvent>(orbit.begin(), orbit.end());
  }
  return event_type_data;
}

}  // namespace

/// \brief Construct SemiGrandCanonicalEventData
SemiGrandCanonicalEventData::SemiGrandCanonicalEventData(
    std::shared_ptr<system_type> system, state_type const &state,
    monte::OccLocation const &occ_location,
    std::vector<monte::OccSwap> const &grand_canonical_swaps,
    std::shared_ptr<SemiGrandCanonicalPotential> potential) {
  // Make OccEvents from SemiGrandCanonical swaps
  // key: event_type_name, value: symmetrically equivalent events
  std::map<std::string, OccEventTypeData> event_type_data =
      _make_event_type_data(system, state, grand_canonical_swaps);

  prim_event_list = clexmonte::make_prim_event_list(event_type_data);

  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, {"formation_energy"});

  // TODO: rejection-kmc option does not require impact table
  event_list = clexmonte::make_complete_event_list(
      prim_event_list, prim_impact_info_list, occ_location);

  // Construct CompleteEventCalculator
  event_calculator = std::make_shared<CompleteEventCalculator>(
      potential, prim_event_list, event_list.events);
}

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM
