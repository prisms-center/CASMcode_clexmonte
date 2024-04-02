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
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccLocation.hh"
#include "casm/monte/io/json/ValueMap_json_io.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;
typedef std::mt19937_64 engine_type;

monte::ValueMap from_variant_type(
    std::variant<monte::ValueMap, nlohmann::json, py::none> const &x) {
  if (x.index() == 0) {
    return std::get<0>(x);
  } else if (x.index() == 1) {
    jsonParser json{static_cast<const nlohmann::json &>(std::get<1>(x))};
    monte::ValueMap values;
    from_json(values, json);
    return values;
  } else if (x.index() == 2) {
    return monte::ValueMap{};
  } else {
    throw std::runtime_error("Unknown error converting to monte::ValueMap");
  }
}

clexmonte::state_type make_state(
    clexmonte::config_type const &configuration,
    std::variant<monte::ValueMap, nlohmann::json, py::none> const &conditions,
    std::variant<monte::ValueMap, nlohmann::json, py::none> const &properties) {
  return clexmonte::state_type(configuration, from_variant_type(conditions),
                               from_variant_type(properties));
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_state, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
      Cluster expansion Monte Carlo state
      )pbdoc";
  py::module::import("libcasm.clexulator");
  py::module::import("libcasm.composition");
  py::module::import("libcasm.configuration");
  py::module::import("libcasm.monte");
  py::module::import("libcasm.xtal");

  py::class_<clexmonte::state_type>(m, "MonteCarloState",
                                    R"pbdoc(
      Cluster expansion model state for Monte Carlo simulations

      The MonteCarloState class holds:

      - the current configuration
      - the thermodynamic conditions
      - configuration properties, if calculated by a Monte Carlo calculator.

      .. rubric:: Special Methods

      - MonteCarloState may be copied with `copy.copy` or `copy.deepcopy`.


      )pbdoc")
      .def(py::init<>(&make_state),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          configuration : libcasm.configuration.Configuration
              The initial configuration (microstate).
          conditions : Union[libcasm.monte.ValueMap, dict, None] = None
              The thermodynamic conditions, as a ValueMap. The accepted
              keys and types depend on the Monte Carlo calculation method and
              are documented with the
              :func:`~libcasm.clexmonte.make_conditions_from_value_map`
              function. If None provided, an empty
              :class:`~libcasm.monte.ValueMap` is used.
          properties : Union[libcasm.monte.ValueMap, dict, None] = None
              Current properties of the state, if provided by the Monte Carlo
              calculation method. If None provided, an empty
              :class:`~libcasm.monte.ValueMap` is used.
          )pbdoc",
           py::arg("configuration"), py::arg("conditions") = std::nullopt,
           py::arg("properties") = std::nullopt)
      .def_readwrite("configuration", &clexmonte::state_type::configuration,
                     R"pbdoc(
          libcasm.configuration.Configuration: The configuration
          )pbdoc")
      .def_readwrite("conditions", &clexmonte::state_type::conditions,
                     R"pbdoc(
         libcasm.monte.ValueMap: The thermodynamic conditions
         )pbdoc")
      .def_readwrite("properties", &clexmonte::state_type::properties,
                     R"pbdoc(
         libcasm.monte.ValueMap: Properties of the state, if provided by the \
         Monte Carlo calculation method.
         )pbdoc")
      .def("__copy__",
           [](clexmonte::state_type const &self) {
             return clexmonte::state_type(self);
           })
      .def("__deepcopy__", [](clexmonte::state_type const &self,
                              py::dict) { return clexmonte::state_type(self); })
      .def_static(
          "from_dict",
          [](const nlohmann::json &data,
             std::shared_ptr<config::SupercellSet> supercells) {
            jsonParser json{data};
            InputParser<clexmonte::state_type> parser(json, *supercells);
            std::runtime_error error_if_invalid{
                "Error in libcasm.clexmonte.MonteCarloState.from_dict"};
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            return std::move(*parser.value);
          },
          R"pbdoc(
          Construct MonteCarloState from a Python dict

          Notes
          -----
          - For a description of the format, see `MonteCarloState JSON object (TODO)`_

          .. _`MonteCarloState JSON object (TODO)`: https://prisms-center.github.io/CASMcode_docs/formats/casm/clex/Configuration/#configdof-json-object

          Parameters
          ----------
          data : dict
            A :class:`~libcasm.clexmonte.MonteCarloState` as a dict.
          supercells : libcasm.configuration.SupercellSet
              A :class:`~libcasm.configuration.SupercellSet`, which holds shared
              supercells in order to avoid duplicates.

          Returns
          -------
          state : libcasm.clexmonte.MonteCarloState
              The :class:`~libcasm.clexmonte.MonteCarloState` constructed from the dict.


          )pbdoc",
          py::arg("data"), py::arg("supercells"))
      .def(
          "to_dict",
          [](clexmonte::state_type const &self, bool write_prim_basis) {
            jsonParser json;
            to_json(self, json, write_prim_basis);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent MonteCarloState as a Python dict

          Notes
          -----
          - For a description of the format, see `MonteCarloState JSON object (TODO)`_

          .. _`MonteCarloState JSON object (TODO)`: https://prisms-center.github.io/CASMcode_docs/formats/casm/clex/Configuration/#configdof-json-object

          Parameters
          ----------
          write_prim_basis : bool, default=False
              If True, write DoF values using the prim basis. Default (False)
              is to write DoF values in the standard basis.

          Returns
          -------
          data : json
              The `MonteCarloState reference (TODO) <https://prisms-center.github.io/CASMcode_docs/formats/casm/clex/Configuration/>`_ documents the expected format for MonteCarloState."
          )pbdoc",
          py::arg("write_prim_basis") = false);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
