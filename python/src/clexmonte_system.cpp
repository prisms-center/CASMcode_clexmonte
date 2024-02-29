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
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

} // namespace CASMpy

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
      .def(py::init<std::shared_ptr<xtal::BasicStructure const> const &,
                    composition::CompositionConverter const &, Index>(),
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
          "prim_neighbor_list",
          [](clexmonte::System &m) -> clexulator::PrimNeighborListWrapper {
            return clexulator::PrimNeighborListWrapper(m.prim_neighbor_list);
          },
          R"pbdoc(
          libcasm.clexulator.PrimNeighborList: Neighbor list used for cluster \
          expansions.
          )pbdoc")
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
          py::arg("data"), py::arg("search_path") = std::vector<std::string>());

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
