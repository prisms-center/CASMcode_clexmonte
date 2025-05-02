#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// nlohmann::json binding
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "pybind11_json/pybind11_json.hpp"

// clexmonte
#include "casm/clexmonte/run/StateModifyingFunction.hh"
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
typedef clexmonte::default_engine_type engine_type;

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

clexmonte::StateModifyingFunction make_modifying_function(
    std::string name, std::string description,
    std::function<void(clexmonte::state_type &, monte::OccLocation *)>
        function) {
  if (function == nullptr) {
    throw std::runtime_error(
        "Error constructing StateModifyingFunction: function == nullptr");
  }
  return clexmonte::StateModifyingFunction(name, description, function);
}

std::shared_ptr<clexmonte::LocalOrbitCompositionCalculator>
make_local_orbit_composition_calculator(
    std::vector<std::vector<std::vector<clust::IntegralCluster>>> const
        &_orbits,
    std::vector<int> _orbits_to_calculate, bool _combine_orbits,
    std::shared_ptr<clexulator::PrimNeighborListWrapper> wrapper,
    std::shared_ptr<clexulator::SuperNeighborList> _supercell_nlist,
    xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
    composition::CompositionCalculator const &_composition_calculator,
    clexulator::ConfigDoFValues const *_dof_values) {
  if (wrapper == nullptr || wrapper->prim_neighbor_list == nullptr) {
    throw std::runtime_error(
        "Error constructing LocalOrbitCompositionCalculator: "
        "prim_neighbor_list is not initialized");
  }

  // Convert the orbits to the correct type
  std::vector<std::vector<std::set<clust::IntegralCluster>>> converted_orbits;
  for (auto const &orbit : _orbits) {
    std::vector<std::set<clust::IntegralCluster>> converted_orbit;
    for (auto const &cluster_set : orbit) {
      std::set<clust::IntegralCluster> converted_cluster_set;
      for (auto const &cluster : cluster_set) {
        converted_cluster_set.insert(cluster);
      }
      converted_orbit.push_back(converted_cluster_set);
    }
    converted_orbits.push_back(converted_orbit);
  }

  // Convert the orbits to calculate to the correct type
  std::set<int> converted_orbits_to_calculate;
  for (auto const &orbit : _orbits_to_calculate) {
    converted_orbits_to_calculate.insert(orbit);
  }

  return std::make_shared<clexmonte::LocalOrbitCompositionCalculator>(
      converted_orbits, converted_orbits_to_calculate, _combine_orbits,
      wrapper->prim_neighbor_list, _supercell_nlist, _supercell_index_converter,
      _composition_calculator, _dof_values);
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MAKE_OPAQUE(CASM::clexmonte::StateModifyingFunctionMap);

PYBIND11_MODULE(_clexmonte_state, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
      Cluster expansion Monte Carlo state
      )pbdoc";
  py::module::import("libcasm.clexulator");
  py::module::import("libcasm.clusterography");
  py::module::import("libcasm.composition");
  py::module::import("libcasm.configuration");
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.events");
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
      .def(
          "copy",
          [](clexmonte::state_type const &self) {
            return clexmonte::state_type(self);
          },
          "Create a copy of the MonteCarloState.")
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

  py::class_<clexmonte::StateModifyingFunction>(m, "StateModifyingFunction",
                                                R"pbdoc(
        Functions that can modify a :class:`~libcasm.clexmonte.MonteCarloState`.

        .. rubric:: Special Methods

        A call operator allows running the function to modify a MonteCarloState.
        Typical usage involves getting a StateModifyingFunction from a
        MonteCalculator and then using it to modify a function.

        .. code-block:: Python

            # calculator: MonteCalculator
            # state: MonteCarloState
            funcname = "enforce.composition"
            enforce_composition_f = calculator.modifying_functions[funcname]
            enforce_composition_f(state)

        )pbdoc")
      .def(py::init<>(&make_modifying_function),
           R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          name : str
              Name of the function
          description : str
              Description of the function.
          function : function
              A function that modifies a
              :class:`~libcasm.clexmonte.MonteCarloState` with signature:

              .. code-block:: Python

                  def func(
                      state: libcasm.clexmonte.MonteCarloState,
                      occ_location: Optional[libcasm.monte.events.OccLocation],
                  ):
                      # function body

              The `state` parameter is the state to be modified and the
              optional `occ_location` parameter allows support for updating an
              :class:`libcasm.monte.events.OccLocation` occupation location
              tracker if the function modifies the state's configuration.

          )pbdoc",
           py::arg("name"), py::arg("description"), py::arg("function"))
      .def_readwrite("name", &clexmonte::StateModifyingFunction::name,
                     R"pbdoc(
          str : Name of the analysis function.
          )pbdoc")
      .def_readwrite("description",
                     &clexmonte::StateModifyingFunction::description,
                     R"pbdoc(
          str : Description of the function.
          )pbdoc")
      .def_readwrite("function", &clexmonte::StateModifyingFunction::function,
                     R"pbdoc(
          function : The state modifying function.

          A function that modifies a
          :class:`~libcasm.clexmonte.MonteCarloState` with signature:

          .. code-block:: Python

              def func(
                  state: libcasm.clexmonte.MonteCarloState,
                  occ_location: Optional[libcasm.monte.events.OccLocation],
              ):
                  # function body

          The `state` parameter is the state to be modified and the
          optional `occ_location` parameter allows support for updating an
          :class:`libcasm.monte.events.OccLocation` occupation location
          tracker if the function modifies the state's configuration.
          )pbdoc")
      .def(
          "__call__",
          [](clexmonte::StateModifyingFunction const &f,
             clexmonte::state_type &state,
             monte::OccLocation *occ_location) { f(state, occ_location); },
          R"pbdoc(
          Runs the state modifying function

          Equivalent to calling
          :py::attr:`~libcasm.clexmonte.StateModifyingFunction.function`.

          Parameters
          ----------
          state: libcasm.clexmonte.MonteCarloState
              The state to be modified
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Occupation location tracker, that will be updated if the state's
              configuration is modified, if supported.
          )pbdoc",
          py::arg("state"),
          py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr));

  py::bind_map<clexmonte::StateModifyingFunctionMap>(
      m, "StateModifyingFunctionMap",
      R"pbdoc(
    StateModifyingFunctionMaP stores :class:`~libcasm.clexmonte.StateModifyingFunction` by name.

    Notes
    -----
    StateModifyingFunctionMap is a Dict[str, :class:`~libcasm.clexmonte.StateModifyingFunction`]-like object.
    )pbdoc",
      py::module_local(false));

  py::class_<clexmonte::LocalOrbitCompositionCalculator,
             std::shared_ptr<clexmonte::LocalOrbitCompositionCalculator>>(
      m, "LocalOrbitCompositionCalculator", R"pbdoc(
        Calculate the composition on sites in local-cluster orbits.

        LocalOrbitCompositionCalculator can be constructed independently, but
        are most conveniently obtained by specifying them in the system input
        file and accessing them using
        :func:`~libcasm.clexmonte.System.local_orbit_composition_calculators`.

        .. rubric:: Special Methods

        A call operator exists which is equivalent to
        :func:`LocalOrbitCompositionCalculator.value`.

        )pbdoc")
      .def(py::init<>(&make_local_orbit_composition_calculator),
           R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          local_orbits : list[list[list[libcasm.clusterography.Cluster]]]
              The local cluster orbits, where
              `orbits[e][i][j]` is the `j`-th cluster in the `i`-th orbit of
              equivalent clusters around the `e`-th equivalent phenomenal
              cluster. Can be obtained from
              :func:`libcasm.clexmonte.System.local_basis_set_cluster_info`.
          orbits_to_calculate : list[int]
              The indices (`i`, begin at 0) of the orbits to calculate. The
              indices will be sorted and duplicates will be removed.
          combine_orbits : bool
              If True, calculate the composition of the union of sites in
              all the `orbits_to_calculate`. If False, calculate the composition
              for the set of sites in each orbit separately.
          prim_nlist : libcasm.clexulator.PrimNeighborList
              The primitive neighbor list.
          supercell_nlist : libcasm.clexulator.SuperNeighborList
              The supercell neighbor list.
          supercell_index_converter : libcasm.xtal.SiteIndexConverter
              The supercell index converter.
          composition_calculator : libcasm.composition.CompositionCalculator
              The composition calculator.
          dof_values : Optional[libcasm.clexulator.ConfigDoFValues]
              The DoF values of the configuration to calculate.

          )pbdoc",
           py::arg("local_orbits"), py::arg("orbits_to_calculate"),
           py::arg("combine_orbits"), py::arg("prim_nlist"),
           py::arg("supercell_nlist"), py::arg("supercell_index_converter"),
           py::arg("composition_calculator"), py::arg("dof_values") = nullptr)
      .def(
          "value",
          [](clexmonte::LocalOrbitCompositionCalculator &f,
             Index unitcell_index, Index equivalent_index) {
            return f.value(unitcell_index, equivalent_index);
          },
          R"pbdoc(
          Calculate the local orbit composition

          Parameters
          ----------
          unitcell_index : int
              The index of the unit cell the phenomenal cluster is associated
              with.
          equivalent_index : int
              The index of the equivalent phenomenal cluster about which the
              local orbits are constructed.

          Returns
          -------
          composition : numpy.ndarray
              The composition (as number of each component), as columns of
              a matrix. If `combine_orbits` is True, then there will be a
              single column; if False, there will be a column for each requested
              orbit in the order given by `orbits_to_calculate`.
          )pbdoc",
          py::arg("unitcell_index"), py::arg("equivalent_index"))
      .def(
          "__call__",
          [](clexmonte::LocalOrbitCompositionCalculator &f,
             Index unitcell_index, Index equivalent_index) {
            return f.value(unitcell_index, equivalent_index);
          },
          py::arg("unitcell_index"), py::arg("equivalent_index"))
      .def_property_readonly(
          "orbits_to_calculate",
          [](clexmonte::LocalOrbitCompositionCalculator const &f) {
            return std::vector<int>(f.orbits_to_calculate().begin(),
                                    f.orbits_to_calculate().end());
          },
          R"pbdoc(
          The orbits to calculate

          Returns
          -------
          orbits_to_calculate : list[int]
              The indices of the orbits to calculate.
          )pbdoc")
      .def_property_readonly(
          "combine_orbits",
          [](clexmonte::LocalOrbitCompositionCalculator const &f) {
            return f.combine_orbits();
          },
          R"pbdoc(
          Whether orbits are combined

          Returns
          -------
          combine_orbits : bool
              If True, calculate the composition of the union of sites in
              all the `orbits_to_calculate`. If False, calculate the composition
              for the set of sites in each orbit separately.
          )pbdoc")
      .def_property_readonly(
          "local_orbits_sites",
          [](clexmonte::LocalOrbitCompositionCalculator const &f) {
            std::vector<std::vector<std::vector<xtal::UnitCellCoord>>>
                converted_local_orbits_sites;
            for (auto const &equiv_orbits : f.local_orbits_sites()) {
              std::vector<std::vector<xtal::UnitCellCoord>>
                  converted_equiv_orbits;
              for (auto const &orbit_sites : equiv_orbits) {
                std::vector<xtal::UnitCellCoord> converted_orbit_sites;
                for (auto const &site : orbit_sites) {
                  converted_orbit_sites.push_back(site);
                }
                converted_equiv_orbits.push_back(converted_orbit_sites);
              }
              converted_local_orbits_sites.push_back(converted_equiv_orbits);
            }
            return converted_local_orbits_sites;
          },
          R"pbdoc(
          The sets of uniques sites used in the calculation

          Returns
          -------
          local_orbits_sites : list[list[list[libcasm.xtal.IntegralSiteCoordinate]]]
              The unique sites in the `orbits_to_calculate`, where
              `orbits[e][i][j]` is the `j`-th site in the `i-th`
              set of sites, for the `e`-th equivalent phenomenal cluster. If
              `combine_orbits` is True, then all sites are combined into the
              same set of sites and `i` only takes value 0; if False, then
              `i` is an index into `orbits_to_calculate`.
          )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
