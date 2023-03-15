#ifndef CASM_clexmonte_events_event_methods
#define CASM_clexmonte_events_event_methods

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"

// debug
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace monte {
class OccLocation;
}

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

/// \brief Construct linear list of events associated with the origin unit
/// cell
template <typename SystemType>
std::vector<PrimEventData> make_prim_event_list(SystemType const &system);

/// \brief Construct linear list of events associated with the origin unit cell
std::vector<PrimEventData> make_prim_event_list(
    std::map<std::string, OccEventTypeData> const &event_type_data);

/// \brief Sets a monte::OccEvent consistent with the PrimEventData and
/// OccLocation
monte::OccEvent &set_event(monte::OccEvent &event,
                           PrimEventData const &prim_event_data,
                           xtal::UnitCell const &translation,
                           monte::OccLocation const &occ_location);

// --- Inline definitions ---

/// \brief Make event required update neighborhood
template <typename SystemType>
EventImpactInfo make_event_impact_info(
    SystemType const &system, PrimEventData const &prim_event_data,
    std::vector<std::string> const &clex_names,
    std::vector<std::string> const &multiclex_names) {
  std::cout << "make_event_impact_info: begin" << std::endl;
  OccEventTypeData const &event_type_data =
      get_event_type_data(system, prim_event_data.event_type_name);
  std::cout << "make_event_impact_info: 1" << std::endl;
  clust::IntegralCluster phenom = make_cluster(prim_event_data.event);

  EventImpactInfo impact;
  impact.phenomenal_sites = phenom.elements();

  // add local basis set dependence
  std::cout << "make_event_impact_info: 2" << std::endl;
  if (!event_type_data.local_multiclex_name.empty()) {
    std::cout << "make_event_impact_info: 3" << std::endl;
    LocalMultiClexData const &local_multiclex_data =
        get_local_multiclex_data(system, event_type_data.local_multiclex_name);
    impact.required_update_neighborhood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index);
  }

  // include impact neighborhood to include clex
  std::cout << "make_event_impact_info: 4" << std::endl;
  for (auto const &name : clex_names) {
    ClexData const &clex_data = get_clex_data(system, name);
    expand(phenom, impact.required_update_neighborhood, clex_data.cluster_info,
           clex_data.coefficients);
  }

  // include impact neighborhood to include multiclex
  std::cout << "make_event_impact_info: 5" << std::endl;
  for (auto const &name : multiclex_names) {
    MultiClexData const &multiclex_data = get_multiclex_data(system, name);
    for (auto const &coeffs : multiclex_data.coefficients) {
      expand(phenom, impact.required_update_neighborhood,
             multiclex_data.cluster_info, coeffs);
    }
  }

  std::cout << "make_event_impact_info: end" << std::endl;
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
