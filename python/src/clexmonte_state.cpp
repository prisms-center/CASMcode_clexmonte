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
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexmonte/system/io/json/System_json_io.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSet.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_state, m) {
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

  py::class_<clexmonte::Configuration>(m, "MonteCarloConfiguration",
                                       R"pbdoc(
      Cluster expansion model configuration for Monte Carlo simulations

      The MonteCarloConfiguration class is a slightly simplified version
      of the :class:`libcasm.configuration.Configuration` data structure, which
      does not include supercell symmetry data.

      Note that this class provides direct access to
      :class:`~libcasm.clexulator.ConfigDoFValues`, which must be used
      carefully. For use with cluster expansion and order parameter calculators,
      DoF values should be in the prim basis and the dimensions of its arrays
      must remain consistent with the prim and supercell transformation matrix.

      Default DoF values can be constructed using:

      - :func:`~libcasm.clexulator.make_default_config_dof_values`: For ConfigDoFValues in the prim basis.
      - :func:`~libcasm.clexulator.make_default_standard_config_dof_values`: For ConfigDoFValues in the standard basis.

      Conversions of DoF values between prim and standard bases, if necessary,
      can be done using:

      - :func:`~libcasm.clexulator.from_standard_values`: Copy ConfigDoFValues and convert from the standard basis to the prim basis.
      - :func:`~libcasm.clexulator.to_standard_values`: Copy ConfigDoFValues and convert from the prim basis to the standard basis.


      .. rubric:: Special Methods

      - MonteCarloConfiguration may be copied with `copy.copy` or `copy.deepcopy`.


      )pbdoc")
      .def(py::init<Eigen::Matrix3l const &,
                    clexulator::ConfigDoFValues const &>(),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
              The transformation matrix, T, relating the superstructure lattice
              vectors, S, to the unit structure lattice vectors, L, according to
              ``S = L @ T``, where S and L are shape=(3,3)  matrices with
              lattice vectors as columns.
          dof_values : libcasm.clexulator.ConfigDoFValues
              The initial occupation, local, and global degree of freedom (DoF)
              values. The dimensions of `dof_values` must be consistent with
              `transformation_matrix_to_super` and the prim. For use with
              cluster expansion and order parameter calculators, DoF values
              should be in the prim basis.
          )pbdoc",
           py::arg("transformation_matrix_to_super"), py::arg("dof_values"))
      .def_readwrite("transformation_matrix_to_super",
                     &clexmonte::Configuration::transformation_matrix_to_super,
                     R"pbdoc(
          numpy.ndarray[numpy.int64[3, 3]]: The supercell transformation matrix

          The transformation matrix, T, relating the superstructure lattice
          vectors, S, to the unit structure lattice vectors, L, according to
          ``S = L @ T``, where S and L are shape=(3,3)  matrices with
          lattice vectors as columns.
          )pbdoc")
      .def_readwrite("dof_values", &clexmonte::Configuration::dof_values,
                     R"pbdoc(
          libcasm.clexulator.ConfigDoFValues: The degree of freedom (DoF) \
          values.

          Note that the dimensions of `dof_values` must remain consistent with
          `transformation_matrix_to_super`. For use with cluster expansion and
          order parameter calculators, DoF values should be in the prim basis.
          )pbdoc")
      .def(
          "dof_values_ptr",
          [](clexmonte::Configuration &self) { return &self.dof_values; },
          "Return a pointer to the ConfigDoFValues")
      .def("__copy__",
           [](clexmonte::Configuration const &self) {
             return clexmonte::Configuration(self);
           })
      .def("__deepcopy__",
           [](clexmonte::Configuration const &self, py::dict) {
             return clexmonte::Configuration(self);
           })
      .def_static(
          "from_dict",
          [](const nlohmann::json &data) {
            jsonParser json{data};
            InputParser<clexmonte::Configuration> parser(json);
            std::runtime_error error_if_invalid{
                "Error in libcasm.clexmonte.MonteCarloConfiguration.from_dict"};
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            return std::move(*parser.value);
          },
          R"pbdoc(
          Construct MonteCarloConfiguration from a Python dict

          Notes
          -----
          - This function does not convert ConfigDoFValues between prim and standard bases, it reads values as they are.
          - Conversions, if necessary, must be done after construction using:

              - :func:`~libcasm.clexulator.from_standard_values`: Copy ConfigDoFValues and convert from the standard basis to the prim basis.
              - :func:`~libcasm.clexulator.to_standard_values`: Copy ConfigDoFValues and convert from the prim basis to the standard basis.

          - For a description of the format, see `ConfigDoF JSON object`_

          .. _`ConfigDoF JSON object`: https://prisms-center.github.io/CASMcode_docs/formats/casm/clex/Configuration/#configdof-json-object

          )pbdoc",
          py::arg("data"))
      .def(
          "to_dict",
          [](clexmonte::Configuration const &self) {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent MonteCarloConfiguration as a Python dict

          Notes
          -----
          - This function does not convert ConfigDoFValues between prim and standard bases, it writes values as they are.
          - Conversions, if necessary, must be done beforehand using:

              - :func:`~libcasm.clexulator.from_standard_values`: Copy ConfigDoFValues and convert from the standard basis to the prim basis.
              - :func:`~libcasm.clexulator.to_standard_values`: Copy ConfigDoFValues and convert from the prim basis to the standard basis.

          - For a description of the format, see `ConfigDoF JSON object`_

          .. _`ConfigDoF JSON object`: https://prisms-center.github.io/CASMcode_docs/formats/casm/clex/Configuration/#configdof-json-object

          )pbdoc")
      .def_static(
          "from_config_with_sym_info",
          [](config::Configuration const &config) -> clexmonte::Configuration {
            return clexmonte::Configuration(
                config.supercell->superlattice.transformation_matrix_to_super(),
                config.dof_values);
          },
          R"pbdoc(
          Construct from a :class:`libcasm.configuration.Configuration`

          This is equivalent to:

          .. code-block:: Python

              # config: libcasm.configuration.Configuration
              mc_config = clexmonte.Configuration(
                  transformation_matrix_to_super=config.supercell().transformation_matrix_to_super(),
                  dof_values=config.dof_values(),
              )

          Parameters
          ----------
          config : libcasm.configuration.Configuration
              A :class:`libcasm.configuration.Configuration` instance.
          )pbdoc",
          py::arg("config"))
      .def(
          "to_config_with_sym_info",
          [](clexmonte::Configuration const &self,
             std::optional<std::shared_ptr<config::Prim const>> prim,
             std::optional<std::shared_ptr<config::SupercellSet>> supercells)
              -> config::Configuration {
            if (supercells.has_value() && supercells.value()) {
              config::SupercellSet &supercell_set = *supercells.value();
              std::shared_ptr<config::Supercell const> supercell =
                  supercell_set.insert(self.transformation_matrix_to_super)
                      .first->supercell;
              return config::Configuration(supercell, self.dof_values);
            } else if (prim.has_value() && prim.value()) {
              std::shared_ptr<config::Supercell const> supercell =
                  std::make_shared<config::Supercell>(
                      prim.value(), self.transformation_matrix_to_super);
              return config::Configuration(supercell, self.dof_values);
            } else {
              throw std::runtime_error(
                  "Error in "
                  "libcasm.clexmonte.MonteCarloConfiguration.to_config_with_"
                  "sym_info: "
                  "One of `prim` or `supercells` must be provided.");
            }
          },
          R"pbdoc(
          Convert to a :class:`libcasm.configuration.Configuration`

          This is equivalent to:

          .. code-block:: Python

              # config: libcasm.configuration.Configuration
              mc_config = clexmonte.Configuration(
                  transformation_matrix_to_super=config.supercell().transformation_matrix_to_super(),
                  dof_values=config.dof_values(),
              )

          Parameters
          ----------
          prim : Optional[libcasm.configuration.Prim] = None
              A :class:`libcasm.configuration.Prim` instance. Only onee of
              `prim` or `supercells` must be provided.
          supercells : Optional[libcasm.configuration.SupercellSet] = None
              A :class:`~libcasm.configuration.SupercellSet`, which holds shared
              supercells in order to avoid duplicates.
          )pbdoc",
          py::arg("prim") = std::nullopt, py::arg("supercells") = std::nullopt);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
