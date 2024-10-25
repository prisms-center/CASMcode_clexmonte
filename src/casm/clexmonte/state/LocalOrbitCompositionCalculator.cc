#include "casm/clexmonte/state/LocalOrbitCompositionCalculator.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/composition/CompositionCalculator.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/monte/sampling/Sampler.hh"

// debug
// #include "casm/clexmonte/definitions.hh"
// #include "casm/clexmonte/state/Configuration.hh"
// #include "casm/clexmonte/system/System.hh"
// #include "casm/configuration/Configuration.hh"
// #include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

/// \brief Calculate the number of possible unique local orbit compositions
///
/// \param comp_calculator The composition calculator
/// \param orbits Stores clusters for each orbit, for each equivalent
/// phenomenal cluster:
///     orbits[equivalent_index][orbit_index] ->
///     std::set<clust::IntegralCluster>
/// \param max_size The maximum number of possible compositions to find
///
/// \return {found_all, found_all ? possible_compositions.size() | max_size}
PossibleLocalOrbitCompositions::PossibleLocalOrbitCompositions(
    composition::CompositionCalculator const &comp_calculator,
    std::vector<std::vector<std::set<clust::IntegralCluster>>> const &orbits,
    Index _max_possible_occupations)
    : max_possible_compositions(_max_possible_occupations),
      found_all(false),
      n_components(comp_calculator.components().size()),
      n_orbits(0) {
  if (orbits.size() == 0) {
    found_all = true;
    return;
  }
  n_orbits = orbits[0].size();

  // Build the counter over possible occupations:
  std::vector<int> initial;
  std::vector<int> final;
  std::vector<int> increment;
  std::vector<int> sublattice_index;
  std::vector<int> orbit_index;

  std::vector<std::vector<Index>> converter =
      composition::make_occ_index_to_component_index_converter(
          comp_calculator.components(), comp_calculator.allowed_occs());

  int i_orbit = 0;
  for (auto const &orbit : orbits[0]) {
    // Get the unique sites in the orbit
    std::set<xtal::UnitCellCoord> sites;
    for (auto const &cluster : orbit) {
      for (auto const &site : cluster) {
        sites.insert(site);
      }
    }

    // Add the sites to the counter
    for (auto const &site : sites) {
      Index b = site.sublattice();
      Index n_occupants = converter[b].size();
      initial.push_back(0);
      final.push_back(n_occupants);
      increment.push_back(1);
      sublattice_index.push_back(b);
      orbit_index.push_back(i_orbit);
    }
    ++i_orbit;
  }

  // Counter over possible occupations on sites in the orbits under
  // consideration
  EigenCounter<Eigen::VectorXi> counter(
      to_VectorXi(initial), to_VectorXi(final), to_VectorXi(increment));

  // Holds the number of each component by orbit for one possible
  // occupation: num_each_component_by_orbit(component_index, orbit_index)
  Eigen::MatrixXi num_each_component_by_orbit =
      Eigen::MatrixXi::Zero(n_components, n_orbits);

  do {
    Eigen::VectorXi const &occ = counter.current();
    num_each_component_by_orbit.setZero();
    for (Index i = 0; i < occ.size(); ++i) {
      num_each_component_by_orbit(converter[sublattice_index[i]][occ[i]],
                                  orbit_index[i])++;
    }

    // Col-major unrolling of the matrix
    Eigen::VectorXl v = num_each_component_by_orbit.reshaped().cast<long>();
    possible.insert(v);

    // If the maximum size of the set of possible occupations has been
    // exceeded, then erase that last one and stop searching
    if (possible.size() > max_possible_compositions) {
      possible.erase(v);
      found_all = false;
      return;
    }

  } while (++counter);

  // All possible occupations have been found
  found_all = true;
}

/// \brief Make a component name, which describes a composition matrix using
/// orbit indices and the component names
///
/// Example output:
/// - If not combining orbits: "[{orbit:1, Va:1, B:1, C:2}, {orbit:2, Va:1, B:2,
/// C:1}]"
/// - If combining orbits: "{orbits:[1,2], Va:2, B:3, C:3}}"
///
/// \param num_each_component_by_orbit The composition matrix being described
/// \param orbits_to_calculate The orbit whose composition is given by each
///     column of the matrix (if not combining orbits), or by the single column
///     (if combining orbits)
/// \param components The names of the components, corresponding to rows of the
///     composition matrix
/// \param combine_orbits If true, the sites in the union of the orbits in
///     `orbits_to_calculate` are considered and the composition is stored in a
///     single column. If false, the sites in each orbit in
///     `orbits_to_calculate` are considered and the composition is stored in a
///     column for each orbit.
/// \return The component name corresponding to the composition matrix
std::string make_component_name(
    Eigen::MatrixXi const &num_each_component_by_orbit,
    std::set<int> orbits_to_calculate,
    std::vector<std::string> const &components, bool combine_orbits) {
  if (combine_orbits) {
    if (num_each_component_by_orbit.cols() != 1) {
      throw std::runtime_error(
          "make_component_name: num_each_component_by_orbit.cols() != 1");
    }

    jsonParser json;
    json["orbits"] = orbits_to_calculate;
    for (Index i = 0; i < num_each_component_by_orbit.rows(); ++i) {
      json[components[i]] = num_each_component_by_orbit(i, 0);
    }
    std::stringstream ss;
    ss << json;
    return ss.str();
  } else {
    if (orbits_to_calculate.size() != num_each_component_by_orbit.cols()) {
      throw std::runtime_error(
          "make_component_name: orbits_to_calculate.size() != "
          "num_each_component_by_orbit.cols()");
    }

    jsonParser json = jsonParser::array();
    int i_col = 0;
    for (int orbit_index : orbits_to_calculate) {
      jsonParser tjson;
      tjson["orbit"] = orbit_index;
      for (Index i = 0; i < num_each_component_by_orbit.rows(); ++i) {
        tjson[components[i]] = num_each_component_by_orbit(i, i_col);
      }
      json.push_back(tjson);
    }
    std::stringstream ss;
    ss << json;
    return ss.str();
  }
}

// LocalOrbitCompositionCalculator::LocalOrbitCompositionCalculator(
//     std::shared_ptr<system_type> _system, std::string _event_type_name,
//     std::set<int> _orbits_to_calculate)
//     : m_system(_system),
//       m_event_type_name(_event_type_name),
//       m_orbits_to_calculate(_orbits_to_calculate) {
//   // Make m_occ_index_to_component_index_converter
//   auto const &composition_calculator = get_composition_calculator(*m_system);
//   m_occ_index_to_component_index_converter =
//       composition::make_occ_index_to_component_index_converter(
//           composition_calculator.components(),
//           composition_calculator.allowed_occs());
//
//   // Setup m_num_each_component_by_orbit and validate orbits_to_calculate
//   auto cluster_info =
//       get_local_basis_set_cluster_info(*m_system, m_event_type_name);
//   int n_orbits = 0;
//   if (cluster_info->orbits.size() > 0) {
//     n_orbits = cluster_info->orbits[0].size();
//   }
//
//   for (int orbit_index : m_orbits_to_calculate) {
//     if (orbit_index < 0 || orbit_index >= n_orbits) {
//       std::stringstream msg;
//       msg << "Error in LocalOrbitCompositionCalculator: "
//           << "orbit_to_calculate=" << orbit_index << " out of range [0,"
//           << n_orbits << ").";
//       throw std::runtime_error(msg.str());
//     }
//   }
//
//   m_num_each_component_by_orbit.resize(
//       composition_calculator.components().size(),
//       m_orbits_to_calculate.size());
//   m_num_each_component_by_orbit.setZero();
// }
//
///// \brief Reset pointer to state currently being calculated
// void LocalOrbitCompositionCalculator::set(state_type const *state) {
//   // supercell-specific
//   m_state = state;
//   if (m_state == nullptr) {
//     throw std::runtime_error(
//         "Error setting LocalOrbitCompositionCalculator state: state is "
//         "empty");
//   }
//
//   // set shell composition calculation data
//   auto cluster_info =
//       get_local_basis_set_cluster_info(*m_system, m_event_type_name);
//
//   // Make m_local_orbits_neighbor_indices:
//   m_local_orbits_neighbor_indices.clear();
//   m_supercell_nlist = get_supercell_neighbor_list(*m_system, *m_state);
//   auto const &convert = get_index_conversions(*m_system, *m_state);
//   auto const &supercell_index_converter = convert.index_converter();
//   for (Index equivalent_index = 0;
//        equivalent_index < cluster_info->orbits.size(); ++equivalent_index) {
//     std::vector<std::set<std::pair<int, int>>> _neighbor_indices_by_orbit;
//     for (auto const &orbit : cluster_info->orbits[equivalent_index]) {
//       std::set<std::pair<int, int>> _neighbor_indices;
//       for (auto const &cluster : orbit) {
//         for (auto const &site : cluster) {
//           Index site_index = supercell_index_converter(site);
//           _neighbor_indices.emplace(
//               m_supercell_nlist->neighbor_index(site_index),
//               site.sublattice());
//         }
//       }
//       _neighbor_indices_by_orbit.emplace_back(std::move(_neighbor_indices));
//     }
//     m_local_orbits_neighbor_indices.emplace_back(
//         std::move(_neighbor_indices_by_orbit));
//   }
// }

/// \brief Constructor - for a single supercell
///
/// \param _orbits The cluster orbits, in order matching a Clexulator, by
///     equivalent index:
///     - The cluster `orbits[equivalent_index][orbit_index][j]` is `j`-th
///       cluster equivalent to the prototype cluster
///       `orbits[equivalent_index][orbit_index][0]` around the
///       `equivalent_index`-th equivalent phenomenal cluster, in the
///       `orbit_index`-th orbit.
/// \param _orbits_to_calculate Orbits to calculate
/// \param _combine_orbits If true, calculate the number of each component for
///      the union of the orbits in `_orbits_to_calculate`. If false, calculate
///      the number of each component for each orbit in `_orbits_to_calculate`
///      individually. If true, the resulting value will be a matrix with a
///      single column, if false, the value will be a matrix with a column for
///      each orbit.
/// \param supercell_nlist Supercell neighbor list
/// \param supercell_index_converter Converter from linear site index to
///     unitcell index and sublattice index
/// \param _dof_values Pointer to the configuration to be calculated (optional).
///     If not provided, the configuration must be set before calling
///     `calculate_num_each_component`.
LocalOrbitCompositionCalculator::LocalOrbitCompositionCalculator(
    std::vector<std::vector<std::set<clust::IntegralCluster>>> const &_orbits,
    std::set<int> _orbits_to_calculate, bool _combine_orbits,
    std::shared_ptr<clexulator::PrimNeighborList> _prim_nlist,
    std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
    xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
    composition::CompositionCalculator const &_composition_calculator,
    clexulator::ConfigDoFValues const *_dof_values)
    : m_orbits_to_calculate(_orbits_to_calculate),
      m_combine_orbits(_combine_orbits),
      m_prim_nlist(_prim_nlist),
      m_supercell_nlist(_supercell_nlist) {
  // Make m_occ_index_to_component_index_converter
  m_occ_index_to_component_index_converter =
      composition::make_occ_index_to_component_index_converter(
          _composition_calculator.components(),
          _composition_calculator.allowed_occs());

  // Validate orbits_to_calculate:
  int n_orbits = 0;
  if (_orbits.size() > 0) {
    n_orbits = _orbits[0].size();
  }

  for (int orbit_index : m_orbits_to_calculate) {
    if (orbit_index < 0 || orbit_index >= n_orbits) {
      std::stringstream msg;
      msg << "Error in LocalOrbitCompositionCalculator: "
          << "orbit_to_calculate=" << orbit_index << " out of range [0,"
          << n_orbits << ").";
      throw std::runtime_error(msg.str());
    }
  }

  // Make:
  // - m_local_orbits_neighbor_indices:
  // - m_local_orbits_sites:
  m_local_orbits_neighbor_indices.clear();
  m_local_orbits_sites.clear();
  for (Index equivalent_index = 0; equivalent_index < _orbits.size();
       ++equivalent_index) {
    std::vector<std::set<std::pair<int, int>>> _neighbor_indices_by_orbit;
    std::vector<std::set<xtal::UnitCellCoord>> _sites_by_orbit;
    int i_orbit = 0;
    for (auto const &orbit : _orbits[equivalent_index]) {
      std::set<std::pair<int, int>> _neighbor_indices;
      std::set<xtal::UnitCellCoord> _sites;
      for (auto const &cluster : orbit) {
        for (auto const &site : cluster) {
          _sites.insert(site);
          _neighbor_indices.emplace(m_prim_nlist->neighbor_index(site),
                                    site.sublattice());
        }
      }
      _neighbor_indices_by_orbit.emplace_back(std::move(_neighbor_indices));
      _sites_by_orbit.emplace_back(std::move(_sites));
      ++i_orbit;
    }
    m_local_orbits_neighbor_indices.emplace_back(
        std::move(_neighbor_indices_by_orbit));
    m_local_orbits_sites.emplace_back(std::move(_sites_by_orbit));
  }

  // Combine orbits if requested
  if (m_combine_orbits) {
    for (int i_equivalent = 0;
         i_equivalent < m_local_orbits_neighbor_indices.size();
         ++i_equivalent) {
      auto &indices = m_local_orbits_neighbor_indices[i_equivalent];
      auto &sites = m_local_orbits_sites[i_equivalent];
      for (int i_orbit = 1; i_orbit < indices.size(); ++i_orbit) {
        indices[0].insert(indices[i_orbit].begin(), indices[i_orbit].end());
        sites[0].insert(sites[i_orbit].begin(), sites[i_orbit].end());
      }
      indices.resize(1);
      sites.resize(1);
    }
    m_unified_orbits_to_calculate = {0};
  } else {
    m_unified_orbits_to_calculate = m_orbits_to_calculate;
  }

  // Setup m_num_each_component_by_orbit
  m_num_each_component_by_orbit.resize(
      _composition_calculator.components().size(),
      m_unified_orbits_to_calculate.size());
  m_num_each_component_by_orbit.setZero();

  // (Optional) Set configuration to be calculated
  if (_dof_values) {
    set(_dof_values);
  }
}

/// \brief Reset pointer to configuration currently being calculated
void LocalOrbitCompositionCalculator::set(
    clexulator::ConfigDoFValues const *dof_values) {
  m_dof_values = dof_values;
  if (m_dof_values == nullptr) {
    throw std::runtime_error(
        "Error setting LocalOrbitCompositionCalculator dof_values: "
        "dof_values is empty");
  }
}

/// \brief Value at particular unit cell and phenomenal cluster
Eigen::MatrixXi const &LocalOrbitCompositionCalculator::value(
    Index unitcell_index, Index equivalent_index) {
  Eigen::VectorXi const &occupation = m_dof_values->occupation;

  std::vector<Index> const &neighbor_index_to_linear_site_index =
      m_supercell_nlist->sites(unitcell_index);

  // indices[orbit_index] = std::set<std::pair<int, int>>
  std::vector<std::set<std::pair<int, int>>> const &indices =
      m_local_orbits_neighbor_indices[equivalent_index];

  m_num_each_component_by_orbit.setZero();
  int col = 0;
  for (int orbit_index : m_unified_orbits_to_calculate) {
    for (auto const &pair : indices[orbit_index]) {
      int neighbor_index = pair.first;
      int sublattice_index = pair.second;
      int site_index = neighbor_index_to_linear_site_index[neighbor_index];
      int occ_index = occupation(site_index);
      int component_index =
          m_occ_index_to_component_index_converter[sublattice_index][occ_index];
      m_num_each_component_by_orbit.col(col)(component_index) += 1;
    }
    ++col;
  }
  return m_num_each_component_by_orbit;
}

}  // namespace clexmonte
}  // namespace CASM
