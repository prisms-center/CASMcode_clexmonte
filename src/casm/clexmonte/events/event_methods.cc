#include "casm/clexmonte/events/event_methods.hh"

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

/// \brief Append events to the prim event list
///
/// \param prim_event_list The prim event list, which will be appended to with
///     all the equivalent events, including the reverse event separately if
///     it is not identical to the forward event (i.e. cyclical events).
/// \param event_type_name The event name
/// \param events The vector of equivalent events in an orbit, which must be in
///     order consistent with equivalents info if the equivalent_index will be
///     used for getting the correct local cluster expansion.
void append_to_prim_event_list(
    std::vector<PrimEventData> &prim_event_list, std::string event_type_name,
    std::vector<occ_events::OccEvent> const &events) {
  Index equivalent_index = 0;
  for (occ_events::OccEvent const &equiv : events) {
    // forward
    PrimEventData data;
    data.event_type_name = event_type_name;
    data.equivalent_index = equivalent_index;
    data.is_forward = true;
    data.prim_event_index = prim_event_list.size();
    data.event = equiv;
    auto clust_occupation = make_cluster_occupation(data.event);
    data.sites = clust_occupation.first.elements();
    data.occ_init = clust_occupation.second[0];
    data.occ_final = clust_occupation.second[1];
    prim_event_list.push_back(data);

    occ_events::OccEvent reverse_equiv = copy_reverse(equiv);
    if (reverse_equiv != equiv) {
      PrimEventData rev_data = data;
      rev_data.is_forward = false;
      rev_data.prim_event_index = prim_event_list.size();
      rev_data.event = reverse_equiv;
      rev_data.occ_init = data.occ_final;
      rev_data.occ_final = data.occ_init;
      prim_event_list.push_back(rev_data);
    }
    ++equivalent_index;
  }
}

/// \brief Construct linear list of events associated with the origin unit cell
std::vector<PrimEventData> make_prim_event_list(
    std::map<std::string, OccEventTypeData> const &event_type_data) {
  std::vector<PrimEventData> prim_event_list;
  for (auto const &pair : event_type_data) {
    append_to_prim_event_list(prim_event_list, pair.first, pair.second.events);
  }
  return prim_event_list;
}

/// \brief Sets `linear_site_index` given a `unitcell_index` and
/// `neighbor_index` list
///
/// This uses the SuperNeighborList to set the linear_site_index for an
/// event in a particular unitcell. It requires that neighbor_index is
/// constructed ahead of time.
void set_event_linear_site_index(
    std::vector<Index> &linear_site_index, Index unitcell_index,
    std::vector<int> neighbor_index,
    clexulator::SuperNeighborList const &supercell_nlist) {
  Index n_sites = neighbor_index.size();
  std::vector<Index> const &neighbor_index_to_linear_site_index =
      supercell_nlist.sites(unitcell_index);

  // set e.linear_site_index --- specify sites being transformed
  linear_site_index.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    linear_site_index[i] =
        neighbor_index_to_linear_site_index[neighbor_index[i]];
  }
}

/// \brief Sets `linear_site_index` given a `unitcell_index` and
/// `neighbor_index` list
///
/// This uses the SuperNeighborList to set the linear_site_index for an
/// event in a particular unitcell. It requires that neighbor_index is
/// constructed ahead of time.
void set_event_occ_transform_and_atom_traj(
    monte::OccEvent &event, PrimEventData const &prim_event_data,
    Index unitcell_index, monte::OccLocation const &occ_location) {
  monte::Conversions const &convert = occ_location.convert();

  Index n_sites = prim_event_data.sites.size();
  event.occ_transform.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    Index l = event.linear_site_index[i];
    Index asym = convert.l_to_asym(l);

    monte::OccTransform &transform = event.occ_transform[i];
    transform.mol_id = occ_location.l_to_mol_id(l);
    transform.l = l;
    transform.asym = asym;
    transform.from_species =
        convert.species_index(asym, prim_event_data.occ_init[i]);
    transform.to_species =
        convert.species_index(asym, prim_event_data.occ_final[i]);
  }

  // set e.atom_traj --- specify atom motion
  xtal::UnitCell const &translation =
      convert.unitcell_index_converter()(unitcell_index);
  Index n_atoms = prim_event_data.event.size();
  event.atom_traj.resize(n_atoms);
  for (Index i = 0; i < n_atoms; ++i) {
    occ_events::OccTrajectory const &occ_traj = prim_event_data.event[i];
    if (occ_traj.position.size() != 2) {
      throw std::runtime_error("Error: KMC event trajectories must be size 2.");
    }
    xtal::UnitCellCoord const &from_site =
        occ_traj.position[0].integral_site_coordinate + translation;
    xtal::UnitCell const &from_unitcell = from_site.unitcell();
    xtal::UnitCellCoord const &to_site =
        occ_traj.position[1].integral_site_coordinate + translation;
    xtal::UnitCell const &to_unitcell = to_site.unitcell();

    monte::AtomTraj &atom_traj = event.atom_traj[i];
    atom_traj.from.l = convert.bijk_to_l(from_site);
    atom_traj.from.mol_id = occ_location.l_to_mol_id(atom_traj.from.l);
    atom_traj.from.mol_comp = occ_traj.position[0].atom_position_index;

    atom_traj.to.l = convert.bijk_to_l(to_site);
    atom_traj.to.mol_id = occ_location.l_to_mol_id(atom_traj.to.l);
    atom_traj.to.mol_comp = occ_traj.position[1].atom_position_index;

    atom_traj.delta_ijk = to_unitcell - from_unitcell;
  }
}

/// \brief Sets `linear_site_index` given a `unitcell_index` and
/// `neighbor_index` list
///
/// This uses the SuperNeighborList to set the linear_site_index for an
/// event in a particular unitcell. It requires that neighbor_index is
/// constructed ahead of time.
monte::OccEvent &set_event(
    monte::OccEvent &event, PrimEventData const &prim_event_data,
    Index unitcell_index, monte::OccLocation const &occ_location,
    std::vector<int> neighbor_index,
    clexulator::SuperNeighborList const &supercell_nlist) {
  // set event.new_occ --- specify new site occupation
  event.new_occ = prim_event_data.occ_final;

  // set event.linear_site_index
  set_event_linear_site_index(event.linear_site_index, unitcell_index,
                              neighbor_index, supercell_nlist);

  // set e.occ_transform --- specify change in occupation variable
  // set e.atom_traj --- specify atom motion
  set_event_occ_transform_and_atom_traj(event, prim_event_data, unitcell_index,
                                        occ_location);

  return event;
}

/// \brief Sets a monte::OccEvent consistent with the PrimEventData and
/// OccLocation
///
/// Notes:
/// - This doesn't need the current occupation state, just unchanging indices
///   into OccLocation, so monte::OccEvent can be set once per supercell and
///   does not need to be updated after an event occurs.
monte::OccEvent &set_event(monte::OccEvent &event,
                           PrimEventData const &prim_event_data,
                           xtal::UnitCell const &translation,
                           monte::OccLocation const &occ_location) {
  for (auto const &traj : prim_event_data.event) {
    for (auto const &pos : traj.position) {
      // if (pos.is_in_resevoir) {
      //   throw std::runtime_error(
      //       "Error: KMC events exchanging with the resevoir is not
      //       allowed.");
      // }
      if (!pos.is_atom) {
        throw std::runtime_error(
            "Error: KMC event trajectories must describe individual atoms.");
      }
    }
  }

  Index n_sites = prim_event_data.sites.size();
  Index n_atoms = prim_event_data.event.size();
  monte::Conversions const &convert = occ_location.convert();
  auto const &unitcellcoord_index_converter = convert.index_converter();

  // set e.new_occ --- specify new site occupation
  event.new_occ = prim_event_data.occ_final;

  // set e.linear_site_index --- specify sites being transformed
  event.linear_site_index.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    event.linear_site_index[i] =
        unitcellcoord_index_converter(prim_event_data.sites[i] + translation);
  }

  // set e.occ_transform --- specify change in occupation variable
  event.occ_transform.resize(n_sites);
  for (Index i = 0; i < n_sites; ++i) {
    Index l = event.linear_site_index[i];
    Index asym = convert.l_to_asym(l);

    monte::OccTransform &transform = event.occ_transform[i];
    transform.mol_id = occ_location.l_to_mol_id(l);
    transform.l = l;
    transform.asym = asym;
    transform.from_species =
        convert.species_index(asym, prim_event_data.occ_init[i]);
    transform.to_species =
        convert.species_index(asym, prim_event_data.occ_final[i]);
  }

  // set e.atom_traj --- specify atom motion
  event.atom_traj.resize(n_atoms);
  for (Index i = 0; i < n_atoms; ++i) {
    occ_events::OccTrajectory const &occ_traj = prim_event_data.event[i];
    if (occ_traj.position.size() != 2) {
      throw std::runtime_error("Error: KMC event trajectories must be size 2.");
    }
    xtal::UnitCellCoord const &from_site =
        occ_traj.position[0].integral_site_coordinate + translation;
    xtal::UnitCell const &from_unitcell = from_site.unitcell();
    xtal::UnitCellCoord const &to_site =
        occ_traj.position[1].integral_site_coordinate + translation;
    xtal::UnitCell const &to_unitcell = to_site.unitcell();

    monte::AtomTraj &atom_traj = event.atom_traj[i];
    atom_traj.from.l = convert.bijk_to_l(from_site);
    atom_traj.from.mol_id = occ_location.l_to_mol_id(atom_traj.from.l);
    atom_traj.from.mol_comp = occ_traj.position[0].atom_position_index;

    atom_traj.to.l = convert.bijk_to_l(to_site);
    atom_traj.to.mol_id = occ_location.l_to_mol_id(atom_traj.to.l);
    atom_traj.to.mol_comp = occ_traj.position[1].atom_position_index;

    atom_traj.delta_ijk = to_unitcell - from_unitcell;
  }

  return event;
}

/// \brief Return true if the event is allowed; false otherwise.
bool event_is_allowed(std::vector<Index> const &linear_site_index,
                      clexulator::ConfigDoFValues const &dof_values,
                      PrimEventData const &prim_event_data) {
  int i = 0;
  for (Index l : linear_site_index) {
    if (dof_values.occupation(l) != prim_event_data.occ_init[i]) {
      return false;
    }
    ++i;
  }
  return true;
}

}  // namespace clexmonte
}  // namespace CASM
