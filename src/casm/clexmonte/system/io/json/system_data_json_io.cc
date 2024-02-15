#include "casm/clexmonte/system/io/json/system_data_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/system_data.hh"
#include "casm/clexulator/Clexulator.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse BasisSetClusterInfo from a bspecs.json / eci.json file
///
/// Notes:
/// - This is valid for periodic, not local-cluster orbits
void parse(
    InputParser<BasisSetClusterInfo> &parser, config::Prim const &prim,
    std::map<std::string, std::shared_ptr<clexulator::Clexulator>> basis_sets) {
  BasisSetClusterInfo curr;

  // "bspecs"/"cluster_specs"/"generating_group"
  std::vector<Index> generating_group_indices;
  parser.require(
      generating_group_indices,
      fs::path("bspecs") / "cluster_specs" / "params" / "generating_group");
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

}  // namespace clexmonte
}  // namespace CASM
