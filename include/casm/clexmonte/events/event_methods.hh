#ifndef CASM_clexmonte_events_event_methods
#define CASM_clexmonte_events_event_methods

#include "casm/clexmonte/events/event_data.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"

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
  OccEventTypeData const &event_type_data =
      get_event_type_data(system, prim_event_data.event_type_name);
  LocalMultiClexData const &local_multiclex_data =
      get_local_multiclex_data(system, event_type_data.local_multiclex_name);

  clust::IntegralCluster phenom = make_cluster(prim_event_data.event);

  EventImpactInfo impact;
  impact.phenomenal_sites = phenom.elements();

  // add local basis set dependence
  impact.required_update_neighborhood = get_required_update_neighborhood(
      system, local_multiclex_data, prim_event_data.equivalent_index);

  // include impact neighborhood to include clex
  for (auto const &name : clex_names) {
    ClexData const &clex_data = get_clex_data(system, name);
    expand(phenom, impact.required_update_neighborhood, clex_data.cluster_info,
           clex_data.coefficients);
  }

  // include impact neighborhood to include multiclex
  for (auto const &name : multiclex_names) {
    MultiClexData const &multiclex_data = get_multiclex_data(system, name);
    for (auto const &coeffs : multiclex_data.coefficients) {
      expand(phenom, impact.required_update_neighborhood,
             multiclex_data.cluster_info, coeffs);
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
  std::vector<PrimEventData> prim_event_list;

  Index event_type_index = 0;
  for (auto const &pair : get_event_type_data(system)) {
    Index equivalent_index = 0;
    for (occ_events::OccEvent const &equiv : pair.second.events) {
      // forward
      PrimEventData data;
      data.event_type_name = pair.first;
      data.equivalent_index = equivalent_index;
      data.is_forward = true;
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
        rev_data.event = reverse_equiv;
        rev_data.occ_init = data.occ_final;
        rev_data.occ_final = data.occ_init;
        prim_event_list.push_back(rev_data);
      }
      ++equivalent_index;
    }
    ++event_type_index;
  }
  return prim_event_list;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
