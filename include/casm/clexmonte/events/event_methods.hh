#ifndef CASM_clexmonte_events_event_methods
#define CASM_clexmonte_events_event_methods

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
class SuperNeighborList;
}  // namespace clexulator

namespace monte {
class OccLocation;
}  // namespace monte

namespace clexmonte {

/// \brief Make event required update neighborhood
template <typename SystemType>
EventImpactInfo make_event_impact_info(
    SystemType const &system, PrimEventData const &prim_event_data,
    std::vector<std::string> const &clex_names = {"formation_energy"},
    std::vector<std::string> const &multiclex_names = {});

/// \brief Construct list of event impact neighborhoods
template <typename SystemType>
std::vector<EventImpactInfo> make_prim_impact_info_list(
    SystemType const &system, std::vector<PrimEventData> const &prim_event_list,
    std::vector<std::string> const &clex_names = {"formation_energy"},
    std::vector<std::string> const &multiclex_names = {});

/// \brief Append events to the prim event list
void append_to_prim_event_list(std::vector<PrimEventData> &prim_event_list,
                               std::string event_type_name,
                               std::vector<occ_events::OccEvent> const &events);

/// \brief Construct linear list of events associated with the origin unit
/// cell
template <typename SystemType>
std::vector<PrimEventData> make_prim_event_list(SystemType const &system);

/// \brief Construct linear list of events associated with the origin unit cell
std::vector<PrimEventData> make_prim_event_list(
    std::map<std::string, OccEventTypeData> const &event_type_data);

/// \brief Sets `linear_site_index` given a `unitcell_index` and
/// `neighbor_index` list
void set_event_linear_site_index(
    std::vector<Index> &linear_site_index, Index unitcell_index,
    std::vector<int> neighbor_index,
    clexulator::SuperNeighborList const &supercell_nlist);

/// \brief Sets `event.occ_transform` and `event.atom_traj`
void set_event_occ_transform_and_atom_traj(
    monte::OccEvent &event, PrimEventData const &prim_event_data,
    Index unitcell_index, monte::OccLocation const &occ_location);

/// \brief Sets a monte::OccEvent
monte::OccEvent &set_event(
    monte::OccEvent &event, PrimEventData const &prim_event_data,
    Index unitcell_index, monte::OccLocation const &occ_location,
    std::vector<int> neighbor_index,
    clexulator::SuperNeighborList const &supercell_nlist);

/// \brief Sets a monte::OccEvent consistent with the PrimEventData and
/// OccLocation
monte::OccEvent &set_event(monte::OccEvent &event,
                           PrimEventData const &prim_event_data,
                           xtal::UnitCell const &translation,
                           monte::OccLocation const &occ_location);

/// \brief Return true if the event is allowed; false otherwise.
bool event_is_allowed(std::vector<Index> const &linear_site_index,
                      clexulator::ConfigDoFValues const &dof_values,
                      PrimEventData const &prim_event_data);

// --- Inline definitions ---

/// \brief Make event required update neighborhood
template <typename SystemType>
EventImpactInfo make_event_impact_info(
    SystemType const &system, PrimEventData const &prim_event_data,
    std::vector<std::string> const &clex_names,
    std::vector<std::string> const &multiclex_names) {
  OccEventTypeData const &event_type_data =
      get_event_type_data(system, prim_event_data.event_type_name);
  clust::IntegralCluster phenom = make_cluster(prim_event_data.event);

  EventImpactInfo impact;
  impact.phenomenal_sites = phenom.elements();

  // add local basis set dependence
  if (!event_type_data.local_multiclex_name.empty()) {
    LocalMultiClexData const &local_multiclex_data =
        get_local_multiclex_data(system, event_type_data.local_multiclex_name);
    impact.required_update_neighborhood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index);
  }

  // include impact neighborhood to include clex
  for (auto const &name : clex_names) {
    ClexData const &clex_data = get_clex_data(system, name);
    if (!clex_data.cluster_info) {
      std::stringstream msg;
      msg << "Error in make_event_impact_info: clex '";
      msg << name << "' does not have cluster_info";
      throw std::runtime_error(msg.str());
    }
    expand(phenom, impact.required_update_neighborhood, *clex_data.cluster_info,
           clex_data.coefficients);
  }

  // include impact neighborhood to include multiclex
  for (auto const &name : multiclex_names) {
    MultiClexData const &multiclex_data = get_multiclex_data(system, name);
    if (!multiclex_data.cluster_info) {
      std::stringstream msg;
      msg << "Error in make_event_impact_info: multiclex '";
      msg << name << "' does not have cluster_info";
      throw std::runtime_error(msg.str());
    }
    for (auto const &coeffs : multiclex_data.coefficients) {
      expand(phenom, impact.required_update_neighborhood,
             *multiclex_data.cluster_info, coeffs);
    }
  }

  return impact;
}

/// \brief Construct list of event impact neighborhoods
template <typename SystemType>
std::vector<EventImpactInfo> make_prim_impact_info_list(
    SystemType const &system, std::vector<PrimEventData> const &prim_event_list,
    std::vector<std::string> const &clex_names,
    std::vector<std::string> const &multiclex_names) {
  std::vector<EventImpactInfo> prim_impact_info_list;
  for (auto const &data : prim_event_list) {
    prim_impact_info_list.push_back(
        make_event_impact_info(system, data, clex_names, multiclex_names));
  }
  return prim_impact_info_list;
}

/// \brief Construct linear list of events associated with the origin unit cell
template <typename SystemType>
std::vector<PrimEventData> make_prim_event_list(SystemType const &system) {
  return make_prim_event_list(get_event_type_data(system));
}

}  // namespace clexmonte
}  // namespace CASM

#endif
