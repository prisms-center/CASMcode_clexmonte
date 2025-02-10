#include "casm/clexmonte/system/io/json/system_data_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/clusterography/io/json/EquivalentsInfo_json_io.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse BasisSetClusterInfo from a basis.json / eci.json file
///
/// Notes:
/// - This is valid for periodic, not local-cluster orbits
void parse(InputParser<BasisSetClusterInfo> &parser, config::Prim const &prim) {
  BasisSetClusterInfo curr;

  // "bspecs"/"cluster_specs"/"generating_group"
  std::vector<Index> generating_group_indices;
  fs::path generating_group_path =
      fs::path("bspecs") / "cluster_specs" / "generating_group";
  if (parser.self.find_at(generating_group_path) != parser.self.end()) {
    parser.require(generating_group_indices, generating_group_path);
  } else {
    fs::path generating_group_path_v1 =
        fs::path("bspecs") / "cluster_specs" / "params" / "generating_group";
    if (parser.self.find_at(generating_group_path_v1) != parser.self.end()) {
      parser.require(generating_group_indices, generating_group_path_v1);
    } else {
      parser.insert_error("generating_group",
                          "A 'generating_group' array is required");
      return;
    }
  }
  if (!parser.valid()) {
    return;
  }

  std::vector<xtal::UnitCellCoordRep> generating_rep;
  for (Index fg_index : generating_group_indices) {
    generating_rep.push_back(
        prim.sym_info.unitcellcoord_symgroup_rep[fg_index]);
  }

  // "orbits"/<i>/"prototype"
  std::vector<clust::IntegralCluster> prototypes;
  if (!parser.self.contains("orbits") || !parser.self["orbits"].is_array()) {
    parser.insert_error("orbits", "An 'orbits' array is required");
    return;
  }
  auto begin = parser.self["orbits"].begin();
  auto end = parser.self["orbits"].end();
  Index orbit_index = 0;
  for (auto it = begin; it != end; ++it) {
    fs::path orbit_path = fs::path("orbits") / std::to_string(orbit_index);

    // "orbits"/<i>/"prototype"
    clust::IntegralCluster prototype;
    parser.require<clust::IntegralCluster>(prototype, orbit_path / "prototype",
                                           *prim.basicstructure);
    prototypes.push_back(prototype);

    // populate function_to_orbit_index
    if (!it->contains("cluster_functions") ||
        !(*it)["cluster_functions"].is_array()) {
      parser.insert_error(orbit_path / "cluster_functions",
                          "A 'cluster_functions' array is required");
      return;
    }
    for (Index j = 0; j < (*it)["cluster_functions"].size(); ++j) {
      curr.function_to_orbit_index.push_back(orbit_index);
    }
    ++orbit_index;
  }

  // generate orbits
  for (auto const &prototype : prototypes) {
    curr.orbits.push_back(make_prim_periodic_orbit(prototype, generating_rep));
  }
  parser.value = std::make_unique<BasisSetClusterInfo>(curr);
}

/// \brief Output minimal "equivalents info" to JSON
jsonParser &to_json(EquivalentsInfo const &equivalents_info, jsonParser &json,
                    xtal::BasicStructure const &prim) {
  // Using CASM::clust::EquivalentsInfo here to avoid duplicating the
  // to_json function
  clust::EquivalentsInfo _info;
  _info.phenomenal_clusters = equivalents_info.phenomenal_clusters;
  _info.equivalent_generating_op_indices =
      equivalents_info.equivalent_generating_op_indices;
  to_json(_info, json, prim);
  return json;
}

/// \brief Parse equivalents_info.json
///
/// TODO: document format (see
/// tests/unit/clexmonte/data/clexmonte/system_template.json)
void parse(InputParser<EquivalentsInfo> &parser, config::Prim const &prim) {
  xtal::BasicStructure const &basicstructure = *prim.basicstructure;

  std::vector<Index> equivalent_generating_op_indices;
  parser.require(equivalent_generating_op_indices, "equivalent_generating_ops");

  std::vector<clust::IntegralCluster> phenomenal_clusters;
  if (parser.self.contains("equivalents")) {
    auto begin = parser.self["equivalents"].begin();
    auto end = parser.self["equivalents"].end();
    int i = 0;
    for (auto it = begin; it != end; ++it) {
      auto subparser = parser.subparse<clust::IntegralCluster>(
          fs::path("equivalents") / std::to_string(i) / "phenomenal",
          basicstructure);
      if (subparser->valid()) {
        phenomenal_clusters.push_back(*subparser->value);
      }
      ++i;
    }
  }

  if (equivalent_generating_op_indices.size() != phenomenal_clusters.size()) {
    parser.insert_error("equivalent_generating_ops",
                        "Size mismatch with 'equivalents'");
  }

  if (equivalent_generating_op_indices.size() == 0) {
    parser.insert_error("equivalent_generating_ops", "Size==0");
  }

  if (!parser.valid()) {
    return;
  }

  parser.value = std::make_unique<EquivalentsInfo>(
      prim, phenomenal_clusters, equivalent_generating_op_indices);
}

/// \brief Parse LocalBasisSetClusterInfo from a basis.json / eci.json file
///
/// Notes:
/// - This is valid for local-cluster orbits, not periodic
void parse(InputParser<LocalBasisSetClusterInfo> &parser,
           config::Prim const &prim, EquivalentsInfo const &info) {
  auto curr = std::make_unique<LocalBasisSetClusterInfo>();

  // "bspecs"/"cluster_specs"/"generating_group"
  std::vector<Index> generating_group_indices;
  fs::path generating_group_path =
      fs::path("bspecs") / "cluster_specs" / "generating_group";
  if (parser.self.find_at(generating_group_path) != parser.self.end()) {
    parser.require(generating_group_indices, generating_group_path);
  } else {
    fs::path generating_group_path_v1 =
        fs::path("bspecs") / "cluster_specs" / "params" / "generating_group";
    if (parser.self.find_at(generating_group_path_v1) != parser.self.end()) {
      parser.require(generating_group_indices, generating_group_path_v1);
    } else {
      parser.insert_error("generating_group",
                          "A 'generating_group' array is required");
      return;
    }
  }
  if (!parser.valid()) {
    return;
  }

  // make local cluster generating group rep
  std::vector<xtal::UnitCellCoordRep> generating_rep;
  if (info.phenomenal_clusters.size() > 0) {
    std::vector<xtal::SymOp> cluster_group_elements;
    for (Index i : generating_group_indices) {
      cluster_group_elements.push_back(clust::make_cluster_group_element(
          info.phenomenal_clusters[0],
          prim.basicstructure->lattice().lat_column_mat(),
          prim.sym_info.factor_group->element[i],
          prim.sym_info.unitcellcoord_symgroup_rep[i]));
    }
    generating_rep = sym_info::make_unitcellcoord_symgroup_rep(
        cluster_group_elements, *prim.basicstructure);
  }

  // "orbits"/<i>/"prototype"
  std::vector<clust::IntegralCluster> prototypes;
  if (!parser.self.contains("orbits") || !parser.self["orbits"].is_array()) {
    parser.insert_error("orbits", "An 'orbits' array is required");
    return;
  }
  auto begin = parser.self["orbits"].begin();
  auto end = parser.self["orbits"].end();
  Index orbit_index = 0;
  for (auto it = begin; it != end; ++it) {
    fs::path orbit_path = fs::path("orbits") / std::to_string(orbit_index);

    // "orbits"/<i>/"prototype"
    clust::IntegralCluster prototype;
    parser.require<clust::IntegralCluster>(prototype, orbit_path / "prototype",
                                           *prim.basicstructure);
    prototypes.push_back(prototype);

    // populate function_to_orbit_index
    if (!it->contains("cluster_functions") ||
        !(*it)["cluster_functions"].is_array()) {
      parser.insert_error(orbit_path / "cluster_functions",
                          "A 'cluster_functions' array is required");
      return;
    }
    for (Index j = 0; j < (*it)["cluster_functions"].size(); ++j) {
      curr->function_to_orbit_index.push_back(orbit_index);
    }
    ++orbit_index;
  }

  // generate orbits about first phenomenal cluster
  std::vector<std::set<clust::IntegralCluster>> orbits;
  if (info.phenomenal_clusters.size() != 0) {
    for (auto const &prototype : prototypes) {
      orbits.push_back(make_local_orbit(prototype, generating_rep));
    }
  }

  // generate equivalent orbits
  for (auto const &equivalent_generating_op : info.equivalent_generating_ops) {
    xtal::UnitCellCoordRep unitcellcoord_rep = xtal::make_unitcellcoord_rep(
        equivalent_generating_op, prim.basicstructure->lattice(),
        xtal::symop_site_map(equivalent_generating_op, *prim.basicstructure));
    std::vector<std::set<clust::IntegralCluster>> equiv_orbits;
    for (auto const &orbit : orbits) {
      std::set<clust::IntegralCluster> equiv_orbit;
      for (auto const &cluster : orbit) {
        equiv_orbit.emplace(copy_apply(unitcellcoord_rep, cluster));
      }
      equiv_orbits.push_back(equiv_orbit);
    }
    curr->orbits.push_back(equiv_orbits);
  }

  parser.value = std::move(curr);
}

/// \brief Output LocalOrbitCompositionCalculatorData as JSON
jsonParser &to_json(LocalOrbitCompositionCalculatorData const &data,
                    jsonParser &json) {
  json.put_obj();
  to_json(data.event_type_name, json["event"]);
  to_json(data.local_basis_set_name, json["local_basis_set"]);
  to_json(data.orbits_to_calculate, json["orbits_to_calculate"]);
  to_json(data.combine_orbits, json["combine_orbits"]);
  to_json(data.max_size, json["max_size"]);
  return json;
}

/// \brief Parse LocalOrbitCompositionCalculatorData from JSON
void parse(InputParser<LocalOrbitCompositionCalculatorData> &parser,
           std::string local_basis_set_name) {
  // parse local-orbit composition calculator data
  LocalOrbitCompositionCalculatorData data;
  data.local_basis_set_name = local_basis_set_name;
  data.event_type_name.clear();
  parser.optional(data.event_type_name, "event");
  parser.require(data.orbits_to_calculate, "orbits_to_calculate");
  parser.require(data.combine_orbits, "combine_orbits");
  data.max_size = 10000;
  parser.optional(data.max_size, "max_size");

  if (!parser.valid()) {
    return;
  }

  parser.value =
      std::make_unique<LocalOrbitCompositionCalculatorData>(std::move(data));
}

}  // namespace clexmonte
}  // namespace CASM
