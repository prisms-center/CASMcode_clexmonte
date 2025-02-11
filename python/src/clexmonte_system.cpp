#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// nlohmann::json binding
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "pybind11_json/pybind11_json.hpp"

// clexmonte
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

std::shared_ptr<clexmonte::System> make_system(
    std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
    composition::CompositionConverter const &_composition_converter,
    Index _n_dimensions) {
  return std::make_shared<clexmonte::System>(
      _shared_prim, _composition_converter, _n_dimensions);
}

template <typename T>
std::vector<std::string> get_keys(std::map<std::string, T> map) {
  std::vector<std::string> keys;
  for (auto const &pair : map) {
    keys.emplace_back(pair.first);
  }
  return keys;
}

std::string resolve_basis_set_name(clexmonte::System const &self,
                                   std::string key, std::stringstream &msg) {
  std::string basis_set_name;
  if (self.basis_sets.count(key)) {
    basis_set_name = key;
    msg << "found basis_set=" << key;
  } else if (self.clex_data.count(key)) {
    basis_set_name = self.clex_data.at(key).basis_set_name;
    msg << "found clex=" << key;
  } else if (self.multiclex_data.count(key)) {
    basis_set_name = self.multiclex_data.at(key).basis_set_name;
    msg << "found multiclex=" << key;
  } else {
    msg << "key=" << key
        << " could not be resolved as a basis set, cluster expansion, "
           "or multi-cluster expansion";
    throw std::runtime_error(msg.str());
  }

  if (!self.basis_sets.count(basis_set_name)) {
    msg << ", but a basis set with basis_set_name=" << basis_set_name
        << " does not exist in the system.";
    throw std::runtime_error(msg.str());
  }

  return basis_set_name;
}

std::string resolve_local_basis_set_name(clexmonte::System const &self,
                                         std::string key,
                                         std::stringstream &msg) {
  std::string local_basis_set_name;
  if (self.local_basis_sets.count(key)) {
    local_basis_set_name = key;
    msg << "found local_basis_set=" << key;
  } else if (self.local_clex_data.count(key)) {
    local_basis_set_name = self.local_clex_data.at(key).local_basis_set_name;
    msg << "found local_clex=" << key;
  } else if (self.local_multiclex_data.count(key)) {
    local_basis_set_name =
        self.local_multiclex_data.at(key).local_basis_set_name;
    msg << "found local_multiclex=" << key;
  } else if (self.event_type_data.count(key)) {
    std::string local_multiclex_name =
        self.event_type_data.at(key).local_multiclex_name;
    msg << "found event_type_name=" << key;
    if (self.local_multiclex_data.count(local_multiclex_name)) {
      local_basis_set_name = self.local_multiclex_data.at(local_multiclex_name)
                                 .local_basis_set_name;
      msg << ", and local_multiclex_name=" << local_multiclex_name;
    } else {
      msg << ", but no local_multiclex_name=" << local_multiclex_name;
      throw std::runtime_error(msg.str());
    }
  } else {
    msg << "key=" << key
        << " could not be resolved as a local basis set, local cluster "
           "expansion, local multi-cluster expansion, or event type name";
    throw std::runtime_error(msg.str());
  }

  if (!self.local_basis_sets.count(local_basis_set_name)) {
    msg << ", but a local basis set with local_basis_set_name="
        << local_basis_set_name << " does not exist in the system.";
    throw std::runtime_error(msg.str());
  }

  return local_basis_set_name;
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_system, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
    Cluster expansion Monte Carlo system
    )pbdoc";
  py::module::import("libcasm.clexmonte._clexmonte_state");
  py::module::import("libcasm.clexulator");
  py::module::import("libcasm.clusterography");
  py::module::import("libcasm.composition");
  py::module::import("libcasm.configuration");
  py::module::import("libcasm.monte.events");
  py::module::import("libcasm.occ_events");
  py::module::import("libcasm.xtal");

  py::class_<clexmonte::System, std::shared_ptr<clexmonte::System>>(m, "System",
                                                                    R"pbdoc(
      Cluster expansion model system data

      The System class:

      - stores property calculators,
      - handles input of data that is used by property calculators, such as
        parametric composition axes, order parameter definitions, neighbor
        lists, and cluster expansion basis sets and coefficients.

      )pbdoc")
      .def(py::init<>(&make_system),
           R"pbdoc(
         .. rubric:: Constructor

         Parameters
         ----------
         xtal_prim : libcasm.xtal.Prim
             A :class:`~libcasm.xtal.Prim`
         composition_converter : libcasm.composition.CompositionConverter
             A :class:`~libcasm.composition.CompositionConverter` instance.
         n_dimensions : int = 3
             Dimensionality used for kinetic coefficients.
         )pbdoc",
           py::arg("xtal_prim"), py::arg("composition_converter"),
           py::arg("n_dimensions") = 3)
      .def_property_readonly(
          "xtal_prim",
          [](clexmonte::System const &m)
              -> std::shared_ptr<xtal::BasicStructure const> {
            return m.prim->basicstructure;
          },
          R"pbdoc(
          libcasm.xtal.Prim: Primitive crystal structure and allowed degrees \
          of freedom (DoF).
          )pbdoc")
      .def_readonly("prim", &clexmonte::System::prim,
                    R"pbdoc(
          libcasm.configuration.Prim: Prim with symmetry information.
          )pbdoc")
      .def_readonly("n_dimensions", &clexmonte::System::n_dimensions,
                    R"pbdoc(
          int: Dimensionality used for kinetic coefficients.
          )pbdoc")
      .def_readonly("composition_converter",
                    &clexmonte::System::composition_converter,
                    R"pbdoc(
          libcasm.composition.CompositionConverter: Converter between number of \
          species per unit cell and parametric composition.
          )pbdoc")
      .def_readonly("composition_calculator",
                    &clexmonte::System::composition_calculator,
                    R"pbdoc(
          libcasm.composition.CompositionCalculator: Calculator for total and \
          sublattice compositions from an integer occupation array.
          )pbdoc")
      .def_property_readonly(
          "species_list",
          [](clexmonte::System &m) -> std::vector<xtal::Molecule> const & {
            return m.convert.species_list();
          },
          R"pbdoc(
          list[libcasm.xtal.Occupant]: List of species (including all \
          orientations) allowed in the system.
          )pbdoc")
      .def_property_readonly(
          "prim_neighbor_list",
          [](clexmonte::System &m) -> clexulator::PrimNeighborListWrapper {
            return clexulator::PrimNeighborListWrapper(m.prim_neighbor_list);
          },
          R"pbdoc(
          libcasm.clexulator.PrimNeighborList: Neighbor list used for cluster \
          expansions.
          )pbdoc")
      .def_property_readonly(
          "basis_set_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.basis_sets);
          },
          R"pbdoc(
          Get a list of basis set keys

          Returns
          -------
          keys : list[str]
              A list of basis set keys.
          )pbdoc")
      .def(
          "is_basis_set",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_basis_set(m, key);
          },
          R"pbdoc(
          Check if a basis set calculator exists

          Parameters
          ----------
          key : str
              Basis set name

          Returns
          -------
          clexulator : libcasm.clexulator.Clexulator
              True if basis set calculator exists for `key`.
          )pbdoc",
          py::arg("key"))
      .def_property_readonly(
          "local_basis_set_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.local_basis_sets);
          },
          R"pbdoc(
          Get a list of local basis set keys

          Returns
          -------
          keys : list[str]
              A list of local basis set keys.
          )pbdoc")
      .def(
          "is_local_basis_set",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_local_basis_set(m, key);
          },
          R"pbdoc(
          Check if a local basis set calculator exists

          Parameters
          ----------
          key : str
              Local basis set name

          Returns
          -------
          exists : bool
              True if local basis set calculator exists for `key`.
          )pbdoc",
          py::arg("key"))
      .def(
          "basis_set",
          [](clexmonte::System &m,
             std::string key) -> std::shared_ptr<clexulator::Clexulator> {
            return clexmonte::get_basis_set(m, key);
          },
          R"pbdoc(
          Get a basis set (Clexulator)

          Parameters
          ----------
          key : str
              Basis set name

          Returns
          -------
          clexulator : libcasm.clexulator.LocalClexulator
              The  cluster expansion basis set calculator.
          )pbdoc",
          py::arg("key"))
      .def(
          "basis_set_cluster_info",
          [](clexmonte::System &m, std::string key) -> py::tuple {
            // print messages to sys.stdout, sys.stderr
            py::scoped_ostream_redirect redirect;
            py::scoped_estream_redirect err_redirect;

            auto it = m.basis_set_cluster_info.find(key);
            if (it == m.basis_set_cluster_info.end()) {
              throw std::runtime_error("Basis set cluster info not found: " +
                                       key);
            }
            auto const &cluster_info = it->second;
            std::vector<std::vector<clust::IntegralCluster>> orbits;
            for (auto const &orbit : cluster_info->orbits) {
              orbits.emplace_back(orbit.begin(), orbit.end());
            }
            return py::make_tuple(orbits,
                                  cluster_info->function_to_orbit_index);
          },
          R"pbdoc(
          Get basis set cluster info

          Parameters
          ----------
          key : str
              Basis set name

          Returns
          -------
          orbits : list[list[libcasm.clusterography.Cluster]]
              The orbits of clusters used to generate the basis set, where
              `orbits[i][j]` is the `j`-th cluster in the `i`-th orbit of
              equivalent clusters. Returned as a copy.

          function_to_orbit_index : list[int]
              The value `i = function_to_orbit_index[j]` specifies that the
              `j`-th function in the basis set involves DoF on the `i`-th orbit
              of clusters. Returned as a copy.

          )pbdoc",
          py::arg("key"))
      .def(
          "local_basis_set",
          [](clexmonte::System &m, std::string key)
              -> std::shared_ptr<clexulator::LocalClexulatorWrapper> {
            return std::make_shared<clexulator::LocalClexulatorWrapper>(
                clexmonte::get_local_basis_set(m, key));
          },
          R"pbdoc(
          Get a local basis set (LocalClexulator)

          Parameters
          ----------
          key : str
              Local basis set name

          Returns
          -------
          local_clexulator : libcasm.clexulator.LocalClexulator
              The local cluster expansion basis set calculator.
          )pbdoc",
          py::arg("key"))
      .def(
          "local_basis_set_cluster_info",
          [](clexmonte::System &m, std::string key) -> py::tuple {
            // print messages to sys.stdout, sys.stderr
            py::scoped_ostream_redirect redirect;
            py::scoped_estream_redirect err_redirect;
            auto it = m.local_basis_set_cluster_info.find(key);
            if (it == m.local_basis_set_cluster_info.end()) {
              throw std::runtime_error(
                  "Local basis set cluster info not found: " + key);
            }
            auto const &cluster_info = it->second;
            std::vector<std::vector<std::vector<clust::IntegralCluster>>>
                orbits;
            for (auto const &equiv_orbits : cluster_info->orbits) {
              std::vector<std::vector<clust::IntegralCluster>> _equiv_orbits;
              for (auto const &orbit : equiv_orbits) {
                _equiv_orbits.emplace_back(orbit.begin(), orbit.end());
              }
              orbits.emplace_back(std::move(_equiv_orbits));
            }
            return py::make_tuple(orbits,
                                  cluster_info->function_to_orbit_index);
          },
          R"pbdoc(
          Get local basis set cluster info

          Parameters
          ----------
          key : str
              Local basis set name

          Returns
          -------
          local_orbits : list[list[list[libcasm.clusterography.IntegralCluster]]]
              The orbits of local-clusters used to generate the basis set, where
              `orbits[e][i][j]` is the `j`-th cluster in the `i`-th orbit of
              equivalent clusters around the `e`-th equivalent phenomenal
              cluster. Returned as a copy.

          function_to_orbit_index : list[int]
              The value `i = function_to_orbit_index[j]` specifies that the
              `j`-th function in the basis set involves DoF on the `i`-th orbit
              of clusters. Returned as a copy.

          )pbdoc",
          py::arg("key"))
      //
      .def_property_readonly(
          "clex_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.clex_data);
          },
          R"pbdoc(
          Get a list of cluster expansion keys

          Returns
          -------
          keys : list[str]
              A list of cluster expansion keys.
          )pbdoc")
      .def(
          "is_clex",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_clex_data(m, key);
          },
          R"pbdoc(
          Check if a cluster expansion exists

          Parameters
          ----------
          key : str
              Cluster expansion name

          Returns
          -------
          exists : bool
              True if cluster expansion exists for `key`.
          )pbdoc",
          py::arg("key"))
      .def_property_readonly(
          "multiclex_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.multiclex_data);
          },
          R"pbdoc(
          Get a list of multi-cluster expansion keys

          Returns
          -------
          keys : list[str]
              A list of multi-cluster expansion keys.
          )pbdoc")
      .def(
          "is_multiclex",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_multiclex_data(m, key);
          },
          R"pbdoc(
          Check if a multi-cluster expansion exists

          Parameters
          ----------
          key : str
              Multi-cluster expansion name

          Returns
          -------
          exists : bool
              True if multi-cluster expansion exists for `key`.
          )pbdoc",
          py::arg("key"))
      .def_property_readonly(
          "local_clex_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.local_clex_data);
          },
          R"pbdoc(
          Get a list of local cluster expansion keys

          Returns
          -------
          keys : list[str]
              A list of local cluster expansion keys.
          )pbdoc")
      .def(
          "is_local_clex",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_local_clex_data(m, key);
          },
          R"pbdoc(
          Check if a local cluster expansion exists

          Parameters
          ----------
          key : str
              Local cluster expansion name

          Returns
          -------
          exists : bool
              True if local cluster expansion exists for `key`.
          )pbdoc",
          py::arg("key"))
      .def_property_readonly(
          "local_multiclex_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.local_multiclex_data);
          },
          R"pbdoc(
          Get a list of local multi-cluster expansion keys

          Returns
          -------
          keys : list[str]
              A list of local multi-cluster expansion keys.
          )pbdoc")
      .def(
          "is_local_multiclex",
          [](clexmonte::System &m, std::string key) -> bool {
            return clexmonte::is_local_multiclex_data(m, key);
          },
          R"pbdoc(
          Check if a local multi-cluster expansion exists

          Parameters
          ----------
          key : str
              Local multi-cluster expansion name

          Returns
          -------
          exists : bool
              True if local multi-cluster expansion exists for `key`.
          )pbdoc",
          py::arg("key"))
      //
      .def(
          "clex",
          [](clexmonte::System &m, clexmonte::state_type const &state,
             std::string key) -> std::shared_ptr<clexulator::ClusterExpansion> {
            return clexmonte::get_clex(m, state, key);
          },
          R"pbdoc(
          Get a cluster expansion calculator

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state to be calculated
          key : str
              Cluster expansion name

          Returns
          -------
          clex : libcasm.clexulator.ClusterExpansion
              The cluster expansion calculator for `key`, set to calculate for
              `state`.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      .def(
          "multiclex",
          [](clexmonte::System &m, clexmonte::state_type const &state,
             std::string key)
              -> std::pair<std::shared_ptr<clexulator::MultiClusterExpansion>,
                           std::map<std::string, Index>> {
            return std::make_pair(
                clexmonte::get_multiclex(m, state, key),
                clexmonte::get_multiclex_data(m, key).coefficients_glossary);
          },
          R"pbdoc(
          Get a multi-cluster expansion calculator

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state to be calculated
          key : str
              Multi-cluster expansion name

          Returns
          -------
          multiclex : libcasm.clexulator.MultiClusterExpansion
              The multi-cluster expansion calculator for `key`, set to
              calculate for `state`.
          glossary : dict[str, int]
              The glossary provides the mapping between the property being
              calculated and the index specifying the order in which
              the MultiClusterExpansion stores coefficients and returns
              property values.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      .def(
          "local_clex",
          [](clexmonte::System &m, clexmonte::state_type const &state,
             std::string key)
              -> std::shared_ptr<clexulator::LocalClusterExpansion> {
            return clexmonte::get_local_clex(m, state, key);
          },
          R"pbdoc(
          Get a local cluster expansion

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state to be calculated
          key : str
              Local cluster expansion name

          Returns
          -------
          local_clex : libcasm.clexulator.LocalClusterExpansion
              The local cluster expansion calculator for `key`, set to
              calculate for `state`.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      .def(
          "local_multiclex",
          [](clexmonte::System &m, clexmonte::state_type const &state,
             std::string key)
              -> std::pair<
                  std::shared_ptr<clexulator::MultiLocalClusterExpansion>,
                  std::map<std::string, Index>> {
            return std::make_pair(clexmonte::get_local_multiclex(m, state, key),
                                  clexmonte::get_local_multiclex_data(m, key)
                                      .coefficients_glossary);
          },
          R"pbdoc(
          Get a local multi-cluster expansion

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state to be calculated
          key : str
              Local multi-cluster expansion name

          Returns
          -------
          local_multiclex : libcasm.clexulator.MultiLocalClusterExpansion
              The local multi-cluster expansion calculator for `key`, set to
              calculate for `state`.
          glossary : dict[str, int]
              The glossary provides the mapping between the property being
              calculated and the index specifying the order in which
              the MultiClusterExpansion stores coefficients and returns
              property values.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      .def(
          "basis_set_name",
          [](clexmonte::System const &self, std::string key) {
            std::stringstream msg;
            msg << "Basis set name not found: ";

            return resolve_basis_set_name(self, key, msg);
          },
          R"pbdoc(
          Get the basis set name for a cluster expansion or multi-cluster
          expansion.

          Parameters
          ----------
          key : str
              A cluster expansion name or multi-cluster expansion name.

          Returns
          -------
          basis_set_name : str
              The basis set name for the specified cluster expansion or
              multi-cluster expansion.
          )pbdoc",
          py::arg("key"))
      //
      .def_property_readonly(
          "dof_space_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.dof_spaces);
          },
          R"pbdoc(
          Get a list of DoFSpace / OrderParameter keys

          Returns
          -------
          keys : list[str]
              A list of DoFSpace / OrderParameter keys.
          )pbdoc")
      .def(
          "dof_space",
          [](clexmonte::System &m,
             std::string key) -> std::shared_ptr<clexulator::DoFSpace const> {
            return m.dof_spaces.at(key);
          },
          R"pbdoc(
          Get the DoFSpace for an order parameter calculator

          Parameters
          ----------
          key : str
              The order parameter name

          Returns
          -------
          dof_space : libcasm.clexulator.DoFSpace
              The DoFSpace of the order parameter calculator for `key`.
          )pbdoc",
          py::arg("key"))
      .def(
          "order_parameter",
          [](clexmonte::System &m, clexmonte::state_type const &state,
             std::string key) -> std::shared_ptr<clexulator::OrderParameter> {
            return clexmonte::get_order_parameter(m, state, key);
          },
          R"pbdoc(
          Get an order parameter calculator

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state to be calculated
          key : str
              The order parameter name

          Returns
          -------
          order_parameter : libcasm.clexulator.OrderParameter
              The order parameter calculator for `key`, set to calculate for
              `state`.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      .def_property_readonly(
          "dof_subspace_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.dof_subspaces);
          },
          R"pbdoc(
          Get a list of keys for DoFSpace / OrderParameter with defined
          subspaces

          Returns
          -------
          keys : list[str]
              A list of keys for DoFSpace / OrderParameter with defined
              subspaces
          )pbdoc")
      .def(
          "order_parameter_subspaces",
          [](clexmonte::System &m,
             std::string key) -> std::vector<std::vector<Index>> {
            return m.dof_subspaces.at(key);
          },
          R"pbdoc(
          Get the indices of DoFSpace basis vectors forming subspaces

          Parameters
          ----------
          key : str
              The order parameter name

          Returns
          -------
          order_parameter_subspaces : list[list[int]]
              The array `order_parameter_subspaces[i]` is the indices of the
              DoFSpace basis vectors that form the `i`-th subspace.
          )pbdoc",
          py::arg("key"))
      //
      .def(
          "canonical_swaps",
          [](clexmonte::System &m) -> std::vector<monte::OccSwap> {
            return clexmonte::get_canonical_swaps(m);
          },
          R"pbdoc(
          Get the swap types for canonical Monte Carlo events

          Returns
          -------
          canonical_swaps : list[libcasm.monte.OccSwap]
              The swap types allowed for canonical Monte Carlo events
          )pbdoc")
      .def(
          "semigrand_canonical_swaps",
          [](clexmonte::System &m) -> std::vector<monte::OccSwap> {
            return clexmonte::get_semigrand_canonical_swaps(m);
          },
          R"pbdoc(
          Get the single site swap types for semi-grand canonical Monte Carlo \
          events

          Returns
          -------
          semigrand_canonical_swaps : list[libcasm.monte.OccSwap]
              The single swap types allowed to be proposed for semi-grand
              canonical Monte Carlo events. May be empty.
          )pbdoc")
      .def(
          "semigrand_canonical_multiswaps",
          [](clexmonte::System &m) -> std::vector<monte::MultiOccSwap> {
            return clexmonte::get_semigrand_canonical_multiswaps(m);
          },
          R"pbdoc(
          Get the multi-site swap types for semi-grand canonical Monte Carlo \
          events

          Returns
          -------
          semigrand_canonical_multiswaps : list[libcasm.monte.OccSwap]
              The multi-site swap types for semi-grand canonical Monte Carlo
              events. May be empty.
          )pbdoc")
      // KMC events
      .def_readonly("event_system", &clexmonte::System::event_system, R"pbdoc(
          libcasm.occ_events.OccSystem: Index conversion tables used for KMC events.
          )pbdoc")
      .def(
          "local_basis_set_name",
          [](clexmonte::System const &self, std::string key) {
            std::stringstream msg;
            msg << "Local basis set name not found: ";

            return resolve_local_basis_set_name(self, key, msg);
          },
          R"pbdoc(
          Get the local basis set name for a local cluster expansion,
          local multi-cluster expansion, or event type.

          Parameters
          ----------
          key : str
              A local cluster expansion name, local multi-cluster expansion
              name, or event type name.

          Returns
          -------
          local_basis_set_name : str
              The local basis set name for the specified local cluster
              expansion basis set.
          )pbdoc",
          py::arg("key"))
      .def(
          "equivalents_info",
          [](clexmonte::System const &self, std::string key) {
            std::stringstream msg;
            msg << "Equivalents info not found: ";

            std::string local_basis_set_name =
                resolve_local_basis_set_name(self, key, msg);

            if (!self.equivalents_info.count(local_basis_set_name)) {
              msg << ", but no equivalents_info for local_basis_set_name="
                  << local_basis_set_name;
              throw std::runtime_error(msg.str());
            }
            auto const &info = self.equivalents_info.at(local_basis_set_name);
            return py::make_tuple(info.phenomenal_clusters,
                                  info.equivalent_generating_op_indices,
                                  info.translations);
          },
          R"pbdoc(
          Get the "equivalents_info" for a local cluster expansion basis set

          Parameters
          ----------
          key : str
              The local basis set name, local cluster expansion name, local
              multi-cluster expansion name, or event type name.

          Returns
          -------
          phenomenal_clusters : list[Cluster]
              The phenomenal clusters of the local basis sets

          equivalent_generating_op_indices : list[int]
              Indices of the factor group operations that
              generate the phenomenal clusters from the
              prototype.

          translations: list[np.ndarray]
              The translations, applied after the equivalent generating factor
              group operations, that result in the equivalent phenomenal clusters.

          )pbdoc",
          py::arg("key"))
      .def(
          "occevent_symgroup_rep",
          [](clexmonte::System const &self) {
            return self.occevent_symgroup_rep;
          },
          R"pbdoc(
          Get the group representation for transforming OccEvent

          Returns
          -------
          occevent_symgroup_rep : list[libcasm.occ_events.OccEventRep]
              Group representation for transforming OccEvent.
          )pbdoc")
      .def_property_readonly(
          "event_type_names",
          [](clexmonte::System const &self) {
            return get_keys(self.event_type_data);
          },
          R"pbdoc(
          list[str]: The list of event type names
          )pbdoc")
      .def(
          "prototype_event",
          [](clexmonte::System const &self, std::string event_type_name) {
            return self.event_type_data.at(event_type_name).prototype_event;
          },
          R"pbdoc(
          Get the prototype OccEvent for a particular event type

          Parameters
          ----------
          event_type_name: str
              The name of the event type

          Returns
          -------
          prototype_event : libcasm.occ_events.OccEvent
              The prototype event for the specified event type.
          )pbdoc",
          py::arg("event_type_name"))
      .def(
          "events",
          [](clexmonte::System const &self, std::string event_type_name) {
            return self.event_type_data.at(event_type_name).events;
          },
          R"pbdoc(
          Get a list of the equivalent OccEvent for a particular event type

          Parameters
          ----------
          event_type_name: str
              The name of the event type

          Returns
          -------
          events : list[libcasm.occ_events.OccEvent]
              A list of the equivalent OccEvent for the specified event type. The
              events are ordered consistently with the local cluster expansion given by
              `events_local_multiclex_name`. This means that `events[equivalent_index]`
              is the phenomenal event for the `equivalent_index`-th local cluster basis
              set.
          )pbdoc",
          py::arg("event_type_name"))
      .def(
          "prim_event_list",
          [](clexmonte::System const &self) { return self.prim_event_list; },
          R"pbdoc(
          Get a linear list of all distinct events associated with the origin
          unit cell

          Returns
          -------
          prim_event_list : list[libcasm.clexmonte.PrimEventData]
              A list of the the distinct PrimEventData, include one entry for
              each distinct event associated with the origin unit cell. Includes
              separate entries for symmetrically equivalent events and for
              forward and reverse events if they are distinct.
          )pbdoc")
      // local orbit composition calculators
      .def_property_readonly(
          "local_orbit_composition_calculator_keys",
          [](clexmonte::System &m) -> std::vector<std::string> {
            return get_keys(m.local_orbit_composition_calculator_data);
          },
          R"pbdoc(
          Get a list of keys for local orbit composition calculators

          Returns
          -------
          keys : list[str]
              A list of keys for local orbit composition calculators
          )pbdoc")
      .def(
          "local_orbit_composition_calculator",
          [](clexmonte::System &self, clexmonte::state_type const &state,
             std::string key) {
            auto &system = self;
            auto data = system.local_orbit_composition_calculator_data.at(key);

            auto const &composition_calculator =
                clexmonte::get_composition_calculator(system);
            auto const &orbits = clexmonte::get_local_basis_set_cluster_info(
                                     system, data->local_basis_set_name)
                                     ->orbits;
            auto prim_nlist = system.prim_neighbor_list;
            auto supercell_nlist =
                clexmonte::get_supercell_neighbor_list(system, state);
            auto const &supercell_index_converter =
                clexmonte::get_index_conversions(system, state)
                    .index_converter();
            clexulator::ConfigDoFValues const *dof_values =
                &clexmonte::get_dof_values(state);

            return std::make_shared<clexmonte::LocalOrbitCompositionCalculator>(
                orbits, data->orbits_to_calculate, data->combine_orbits,
                prim_nlist, supercell_nlist, supercell_index_converter,
                composition_calculator, dof_values);
          },
          R"pbdoc(
          Get a local orbit composition calculator

          Standard local orbit composition calculators are generated for each
          KMC event using the local-cluster orbits of the local basis set used
          to parameterize its properties. The standard composition calculators
          are given the following names:

          - "<event_name>-<i>": The composition on the union of sites
            in the `i`-th local-cluster orbit only, where `i` is the linear
            orbit index for orbit being calculated, as a matrix with a
            single column.
          - "<event_name>-all": Calculates the composition for each point
            cluster orbit in separate columns.
          - "<event_name>-all-combined": Calculates the composition for the
            union of sites in all point cluster orbits, as a matrix with a
            single column.

          Custom local orbit composition calculators are generated if they are
          specified using the `local_orbit_composition` attribute for a local
          basis set in the system input file. The `local_orbit_composition`
          is a dict attribute that can be used to specify one or more
          calculators, where the keys are names for the calculators and the
          values have the following format:

          - "orbits_to_calculate": list[int], The indices of the orbits of
            local-clusters to include in the composition calculation. The
            indices (begin at 0) correspond to the linear orbit indices in the
            associated `basis.json` file and the `local_orbits` obtained from
            :func:`System.local_basis_set_cluster_info`.
          - "combine_orbits": bool, If true, the composition will be calculated
            for the union of sites in the orbits specified in
            `orbits_to_calculate`. If false, the composition will be calculated
            for the set of sites in each orbit separately.
          - "event": Optional[str], The name of the type of KMC event that
            this calculator is associated with. If specified, a function named
            `local_orbit_composition.<key>` will be generated as one of the
            standard event data collecting functions that can be used during
            KMC simulations to sample the local orbit composition when this
            type of event is selected. If not specified, the calculator can
            still be accessed, but a selected event data collecting function
            that can be used during KMC simulations will not be constructed.
          - "max_size": Optional[int] = 10000, For selected event data
            collection during KMC simulations, a histogram is constructed for
            the number of occurrences of each encountered local orbit
            composition and this number gives the maximum number of unique
            compositions that are tracked. If additional unique local orbit
            compositions are encountered their counts are added to the
            out-of-range bin.

          Example input:

          .. code-block:: python

              "local_basis_sets": {
                "A_Va_1NN": {
                  "equivalents_info": ...
                  "source": ...
                  "basis": ...
                  "local_orbit_composition": {
                    "A_Va_1NN-1+3": {
                      "event": "A_Va_1NN",
                      "orbits_to_calculate": [1, 3],
                      "combine_orbits": true,
                      "max_size": 10000
                    },
                    "A_Va_1NN-1:3": {
                      "event": "A_Va_1NN",
                      "orbits_to_calculate": [1, 2, 3],
                      "combine_orbits": true,
                      "max_size": 10000
                    }
                  }
                },
                ... other local basis sets ...
              }

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The state in which the local orbit composition is to be calculated
          key : str
              The local orbit composition calculator name.

          Returns
          -------
          local_orbit_composition_calculator : LocalOrbitCompositionCalculator
              The local orbit composition calculator.
          )pbdoc",
          py::arg("state"), py::arg("key"))
      //
      .def_static(
          "from_dict",
          [](const nlohmann::json &data, std::vector<std::string> _search_path,
             bool verbose) {
            // print messages to sys.stdout, sys.stderr
            py::scoped_ostream_redirect redirect;
            py::scoped_estream_redirect err_redirect;

            jsonParser json{data};
            std::vector<fs::path> search_path(_search_path.begin(),
                                              _search_path.end());
            InputParser<clexmonte::System> parser(json, search_path, verbose);
            std::runtime_error error_if_invalid{
                "Error in libcasm.clexmonte.System.from_dict"};
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            std::shared_ptr<clexmonte::System> system(parser.value.release());
            return system;
          },
          R"pbdoc(
          Construct a System from a Python dict.

          Parameters
          ----------
          data: dict
              A Python dict, with a format as specified by the
              `System reference <https://prisms-center.github.io/CASMcode_docs/formats/casm/clexmonte/System/>`_
          search_path: list[str] = []
              Relative file paths included in `data` are searched for relative
              to the paths specified by `search_path`.
          verbose: bool = False
              Print progress statements during parsing for debugging purposes.

          )pbdoc",
          py::arg("data"), py::arg("search_path") = std::vector<std::string>(),
          py::arg("verbose") = false)
      .def("make_default_configuration", &clexmonte::make_default_configuration,
           R"pbdoc(
          Construct a default configuration in a specified supercell

          Parameters
          ----------
          transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
              The transformation matrix, T, relating the superstructure lattice
              vectors, S, to the unit structure lattice vectors, L, according to
              ``S = L @ T``, where S and L are shape=(3,3)  matrices with
              lattice vectors as columns.

          Returns
          -------
          default_configuration : libcasm.configuration.Configuration
              A configuration in the specified supercell, with DoF values
              expressed in the prim basis, initialized to default values (0 for
              occupation indices, 0.0 for all global and local DoF components).
          )pbdoc",
           py::arg("transformation_matrix_to_super"))
      .def("make_default_state", &clexmonte::make_default_state,
           R"pbdoc(
          Construct a default MonteCarloState in a specified supercell

          Parameters
          ----------
          transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
              The transformation matrix, T, relating the superstructure lattice
              vectors, S, to the unit structure lattice vectors, L, according to
              ``S = L @ T``, where S and L are shape=(3,3)  matrices with
              lattice vectors as columns.
          conditions :


          Returns
          -------
          default_state : MonteCarloState
              A state in the specified supercell, with a default configuration
              having DoF values expressed in the prim basis, initialized to
              default values (0 for occupation indices, 0.0 for all global and
              local DoF components), and empty conditions.
          )pbdoc",
           py::arg("transformation_matrix_to_super"))
      .def_readonly("supercells", &clexmonte::System::supercells, R"pbdoc(
          libcasm.configuration.SupercellSet: Shares supercells used for
          MonteCarloState.
          )pbdoc")
      .def_property_readonly(
          "additional_params",
          [](clexmonte::System const &self) {
            return static_cast<nlohmann::json>(self.additional_params);
          },
          R"pbdoc(
          dict: A dict of additional parameters, which may be used to customize
          sampling functions or other purposes.
          )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
