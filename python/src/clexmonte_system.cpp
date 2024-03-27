#include <pybind11/eigen.h>
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

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_system, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Cluster expansion Monte Carlo classes and methods

        libcasm.clexmonte._system
        -------------------------

        The System class:

        - stores property calculators,
        - handle input of data that is used by property calculators, such as
          parametric composition axes, order parameter definitions, neighbor
          lists, and cluster expansion basis sets and coefficients.

    )pbdoc";
  py::module::import("libcasm.clexulator");
  py::module::import("libcasm.composition");
  py::module::import("libcasm.configuration");
  py::module::import("libcasm.monte.events");
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
      //
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
              -> std::shared_ptr<clexulator::MultiClusterExpansion> {
            return clexmonte::get_multiclex(m, state, key);
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
              -> std::shared_ptr<clexulator::MultiLocalClusterExpansion> {
            return clexmonte::get_local_multiclex(m, state, key);
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
          )pbdoc",
          py::arg("state"), py::arg("key"))
      //
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
      //
      .def_static(
          "from_dict",
          [](const nlohmann::json &data,
             std::vector<std::string> _search_path) {
            jsonParser json{data};
            std::vector<fs::path> search_path(_search_path.begin(),
                                              _search_path.end());
            InputParser<clexmonte::System> parser(json, search_path);
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
          )pbdoc",
          py::arg("data"), py::arg("search_path") = std::vector<std::string>())
      .def("make_default_configuration", &clexmonte::make_default_configuration,
           R"pbdoc(
          Construct a default MonteCarloConfiguration in a specified supercell

          Parameters
          ----------
          transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
              The transformation matrix, T, relating the superstructure lattice
              vectors, S, to the unit structure lattice vectors, L, according to
              ``S = L @ T``, where S and L are shape=(3,3)  matrices with
              lattice vectors as columns.

          Returns
          -------
          default_configuration : MonteCarloConfiguration
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

          Returns
          -------
          default_state : MonteCarloState
              A state in the specified supercell, with a default configuration
              having DoF values expressed in the prim basis, initialized to
              default values (0 for occupation indices, 0.0 for all global and
              local DoF components), and empty conditions.
          )pbdoc",
           py::arg("transformation_matrix_to_super"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
