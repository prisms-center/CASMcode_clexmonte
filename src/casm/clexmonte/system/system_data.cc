#include "casm/clexmonte/system/system_data.hh"

#include "casm/configuration/Prim.hh"
#include "casm/configuration/clusterography/impact_neighborhood.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"

namespace CASM {
namespace clexmonte {

/// \brief Expand a required_update_neighborhood based on
///     BasisSetClusterInfo and SparseCoefficients
void expand(clust::IntegralCluster const &phenom,
            std::set<xtal::UnitCellCoord> &required_update_neighborhood,
            BasisSetClusterInfo const &cluster_info,
            clexulator::SparseCoefficients const &coefficients) {
  for (Index function_index : coefficients.index) {
    Index orbit_index = cluster_info.function_to_orbit_index[function_index];
    add_to_flower_neighborhood(phenom, required_update_neighborhood,
                               cluster_info.orbits[orbit_index]);
  }
}

EquivalentsInfo::EquivalentsInfo(
    config::Prim const &_prim,
    std::vector<clust::IntegralCluster> const &_phenomenal_clusters,
    std::vector<Index> const &_equivalent_generating_op_indices)
    : phenomenal_clusters(_phenomenal_clusters),
      equivalent_generating_op_indices(_equivalent_generating_op_indices) {
  if (equivalent_generating_op_indices.size() != phenomenal_clusters.size()) {
    throw std::runtime_error(
        "Error constructing clexmonte::EquivalentsInfo: phenomenal_clusters "
        "and equivalent_generating_op_indices size mismatch");
  }
  if (equivalent_generating_op_indices.size() == 0) {
    throw std::runtime_error(
        "Error constructing clexmonte::EquivalentsInfo: "
        "equivalent_generating_op_indices size==0");
  }

  auto const &fg_element = _prim.sym_info.factor_group->element;
  for (Index i = 0; i < equivalent_generating_op_indices.size(); ++i) {
    Index fg_index = equivalent_generating_op_indices[i];
    if (fg_index < 0 || fg_index >= fg_element.size()) {
      throw std::runtime_error(
          "Error constructing clexmonte::EquivalentsInfo: Invalid "
          "equivalent_generating_op_indices value");
    }

    xtal::SymOp const &factor_group_op =
        _prim.sym_info.factor_group->element[fg_index];
    xtal::UnitCellCoordRep unitcellcoord_rep =
        _prim.sym_info.unitcellcoord_symgroup_rep[fg_index];
    Eigen::Matrix3d const &lat_column_mat =
        _prim.basicstructure->lattice().lat_column_mat();

    // get appropriate translation
    xtal::UnitCell translation = equivalence_map_translation(
        unitcellcoord_rep, phenomenal_clusters[0], phenomenal_clusters[i]);
    translations.push_back(translation);

    // get equivalence map op
    xtal::SymOp translation_op(Eigen::Matrix3d::Identity(),
                               lat_column_mat * translation.cast<double>(),
                               false);
    xtal::SymOp equivalence_map_op = translation_op * factor_group_op;
    equivalent_generating_ops.push_back(equivalence_map_op);
  }
}

/// \brief Make equivalents by applying symmetry.
///
/// Generates equivalents according to:
/// \code
/// Index fg_index = info.equivalent_generating_op_indices[i];
/// equivalents[i] = copy_apply(
///     occevent_symgroup_rep[fg_index], event) + info.translations[i];
/// \endcode
///
/// \param event OccEvent to make equivalents of. Should have the same
///     phenomenal cluster used to generate the local basis set.
/// \param info EquivalentsInfo describing symmetry operations used to generate
///     equivalent local basis sets.
/// \param occevent_symgroup_rep Representation of the prim factor group.
std::vector<occ_events::OccEvent> make_equivalents(
    occ_events::OccEvent const &event, EquivalentsInfo const &info,
    std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep) {
  // generate equivalent events, consistent with the symmetry used to
  // generate the local basis set
  std::vector<occ_events::OccEvent> equivalents;
  for (Index i = 0; i < info.translations.size(); ++i) {
    Index fg_index = info.equivalent_generating_op_indices[i];
    equivalents.push_back(copy_apply(occevent_symgroup_rep[fg_index], event) +
                          info.translations[i]);
  }
  return equivalents;
}

bool is_same_phenomenal_clusters(
    std::vector<occ_events::OccEvent> const &equivalents,
    EquivalentsInfo const &info) {
  for (Index i = 0; i < equivalents.size(); ++i) {
    // check that phenomenal_clusters agree
    clust::IntegralCluster equiv_event_cluster = make_cluster(equivalents[i]);
    equiv_event_cluster.sort();
    clust::IntegralCluster equiv_basis_set_phenom = info.phenomenal_clusters[i];
    equiv_basis_set_phenom.sort();
    if (equiv_event_cluster != equiv_basis_set_phenom) {
      return false;
    }
  }
  return true;
}

}  // namespace clexmonte
}  // namespace CASM
