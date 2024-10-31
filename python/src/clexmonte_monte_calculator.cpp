#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
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

// clexmonte/semigrand_canonical
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/io/json/MonteCalculator_json_io.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

extern "C" {
/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// SemiGrandCanonicalCalculator
CASM::clexmonte::BaseMonteCalculator *make_SemiGrandCanonicalCalculator();

/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// CanonicalCalculator
CASM::clexmonte::BaseMonteCalculator *make_CanonicalCalculator();

/// \brief Returns a clexmonte::BaseMonteCalculator* owning a
/// KineticCalculator
CASM::clexmonte::BaseMonteCalculator *make_KineticCalculator();
}

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// used for libcasm.clexmonte:
typedef std::mt19937_64 engine_type;
typedef monte::RandomNumberGenerator<engine_type> generator_type;
typedef clexmonte::MonteCalculator calculator_type;
typedef clexmonte::MontePotential potential_type;
typedef clexmonte::config_type config_type;
typedef clexmonte::state_type state_type;
typedef clexmonte::statistics_type statistics_type;
typedef clexmonte::System system_type;
typedef monte::SamplingFixture<config_type, statistics_type, engine_type>
    sampling_fixture_type;
typedef clexmonte::sampling_fixture_params_type sampling_fixture_params_type;
typedef clexmonte::run_manager_type<engine_type> run_manager_type;
typedef monte::ResultsAnalysisFunction<config_type, statistics_type>
    analysis_function_type;
typedef monte::ResultsAnalysisFunctionMap<config_type, statistics_type>
    analysis_function_map_type;

clexmonte::MontePotential make_potential(
    std::shared_ptr<clexmonte::MonteCalculator> calculator, state_type &state) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;
  monte::OccLocation *occ_location = nullptr;
  calculator->set_state_and_potential(state, nullptr);
  return calculator->potential();
}

std::shared_ptr<clexmonte::EventDataSummary> make_event_data_summary(
    std::shared_ptr<clexmonte::MonteCalculator> calculator,
    double energy_bin_width, double freq_bin_width, double rate_bin_width) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;

  return std::make_shared<clexmonte::EventDataSummary>(
      calculator->state_data(), calculator->event_data(), energy_bin_width,
      freq_bin_width, rate_bin_width);
}

clexmonte::MonteEventData make_event_data(
    std::shared_ptr<clexmonte::MonteCalculator> calculator, state_type &state,
    std::shared_ptr<engine_type> engine, monte::OccLocation *occ_location) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;
  calculator->set_state_and_potential(state, occ_location);
  if (occ_location == nullptr) {
    calculator->state_data()->owned_occ_location =
        calculator->make_occ_location();
  }
  calculator->set_event_data(engine);
  return calculator->event_data();
}

std::shared_ptr<clexmonte::MonteCalculator>
make_shared_SemiGrandCanonicalCalculator(jsonParser const &params,
                                         std::shared_ptr<system_type> system) {
  std::shared_ptr<RuntimeLibrary> lib = nullptr;
  return clexmonte::make_monte_calculator(
      params, system,
      std::unique_ptr<clexmonte::BaseMonteCalculator>(
          make_SemiGrandCanonicalCalculator()),
      lib);
}

std::shared_ptr<clexmonte::MonteCalculator> make_shared_CanonicalCalculator(
    jsonParser const &params, std::shared_ptr<system_type> system) {
  std::shared_ptr<RuntimeLibrary> lib = nullptr;
  return clexmonte::make_monte_calculator(
      params, system,
      std::unique_ptr<clexmonte::BaseMonteCalculator>(
          make_CanonicalCalculator()),
      lib);
}

std::shared_ptr<clexmonte::MonteCalculator> make_shared_KineticCalculator(
    jsonParser const &params, std::shared_ptr<system_type> system) {
  std::shared_ptr<RuntimeLibrary> lib = nullptr;
  return clexmonte::make_monte_calculator(
      params, system,
      std::unique_ptr<clexmonte::BaseMonteCalculator>(make_KineticCalculator()),
      lib);
}

std::shared_ptr<clexmonte::StateData> make_state_data(
    std::shared_ptr<system_type> system, state_type &state,
    monte::OccLocation *occ_location) {
  return std::make_shared<clexmonte::StateData>(system, &state, occ_location);
}

std::shared_ptr<clexmonte::MonteCalculator> make_monte_calculator(
    std::string method, std::shared_ptr<system_type> system,
    std::optional<nlohmann::json> params) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;
  jsonParser _params = jsonParser::object();
  if (params.has_value()) {
    jsonParser json{static_cast<nlohmann::json const &>(params.value())};
    _params = json;
  }

  if (method == "semigrand_canonical") {
    return make_shared_SemiGrandCanonicalCalculator(_params, system);
  } else if (method == "canonical") {
    return make_shared_CanonicalCalculator(_params, system);
  } else if (method == "kinetic") {
    return make_shared_KineticCalculator(_params, system);
  } else {
    std::stringstream msg;
    msg << "Error in make_monte_calculator: method='" << method
        << "' is not recognized";
    throw std::runtime_error(msg.str());
  }
}

std::shared_ptr<clexmonte::MonteCalculator> make_custom_monte_calculator(
    std::shared_ptr<system_type> system, std::string source,
    std::optional<nlohmann::json> params,
    std::optional<std::string> compile_options,
    std::optional<std::string> so_options,
    std::optional<std::vector<std::string>> search_path) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;
  // fs::path dirpath, std::string calculator_name

  jsonParser _params = jsonParser::object();
  if (params.has_value()) {
    jsonParser json{static_cast<nlohmann::json const &>(params.value())};
    _params = json;
  }

  // Use JSON parser to avoid duplication and give nice error messages
  jsonParser json;
  json["source"] = source;
  if (compile_options.has_value()) {
    json["compile_options"] = compile_options.value();
  }
  if (so_options.has_value()) {
    json["so_options"] = compile_options.value();
  }

  std::vector<fs::path> _search_path;
  if (search_path.has_value()) {
    for (auto const &tpath : search_path.value()) {
      _search_path.emplace_back(tpath);
    }
  }
  InputParser<std::shared_ptr<clexmonte::MonteCalculator>> parser(
      json, system, _params, _search_path);
  std::runtime_error error_if_invalid{
      "Error in libcasm.clexmonte.make_monte_calculator"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  return *parser.value;
}

std::shared_ptr<run_manager_type> monte_calculator_run(
    calculator_type &self, state_type &state,
    std::shared_ptr<run_manager_type> run_manager,
    monte::OccLocation *occ_location) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;

  if (run_manager == nullptr) {
    throw std::runtime_error(
        "Error in MonteCalculator.run: run_manager is None");
  }

  if (run_manager->sampling_fixtures.size() == 0) {
    throw std::runtime_error(
        "Error in MonteCalculator.run: "
        "len(run_manager.sampling_fixtures) == 0");
  }

  // Need to check for an OccLocation
  std::unique_ptr<monte::OccLocation> tmp;
  make_temporary_if_necessary(state, occ_location, tmp, self);

  // run
  self.run(state, *occ_location, run_manager);
  return run_manager;
}

std::shared_ptr<sampling_fixture_type> monte_calculator_run_fixture(
    calculator_type &self, state_type &state,
    sampling_fixture_params_type &sampling_fixture_params,
    std::shared_ptr<engine_type> engine, monte::OccLocation *occ_location) {
  // print errors and warnings to sys.stdout
  py::scoped_ostream_redirect redirect;
  if (!engine) {
    engine = std::make_shared<engine_type>();
    std::random_device device;
    engine->seed(device());
  }
  std::vector<sampling_fixture_params_type> _sampling_fixture_params;
  _sampling_fixture_params.push_back(sampling_fixture_params);
  bool global_cutoff = true;
  std::shared_ptr<run_manager_type> run_manager =
      std::make_shared<run_manager_type>(engine, _sampling_fixture_params,
                                         global_cutoff);
  // run
  monte_calculator_run(self, state, run_manager, occ_location);
  return run_manager->sampling_fixtures.at(0);
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MAKE_OPAQUE(CASM::monte::SamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonSamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::StateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonStateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASMpy::analysis_function_map_type);
PYBIND11_MAKE_OPAQUE(CASM::clexmonte::StateModifyingFunctionMap);
PYBIND11_MAKE_OPAQUE(std::vector<CASM::clexmonte::PrimEventData>);

PYBIND11_MODULE(_clexmonte_monte_calculator, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
    Cluster expansion Monte Carlo implementations
    )pbdoc";
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.events");
  py::module::import("libcasm.monte.sampling");
  py::module::import("libcasm.clexmonte._clexmonte_system");
  py::module::import("libcasm.clexmonte._clexmonte_state");
  py::module::import("libcasm.clexmonte._clexmonte_run_management");

  py::class_<clexmonte::StateData, std::shared_ptr<clexmonte::StateData>>(
      m, "StateData",
      R"pbdoc(
      Access state-specific data used in a Monte Carlo method

      )pbdoc")
      .def(py::init<>(&make_state_data),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        system : libcasm.clexmonte.System
            Cluster expansion model system data.
        state : libcasm.clexmonte.MonteCarloState
            The input state.
        occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state` and it is used and updated during the run.
              If None, no occupant location list is stored. The occupant
              location list is not required for evaluating the potential.
        )pbdoc",
           py::arg("system"), py::arg("state"),
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
      .def_readonly("system", &clexmonte::StateData::system, R"pbdoc(
          System : System data.
          )pbdoc")
      .def_readonly("state", &clexmonte::StateData::state, R"pbdoc(
          Optional[MonteCarloState] : The current state.
          )pbdoc")
      .def_readonly("transformation_matrix_to_super",
                    &clexmonte::StateData::transformation_matrix_to_super,
                    R"pbdoc(
          np.ndarray[np.int64] : The current state's supercell transformation \
          matrix.
          )pbdoc")
      .def_readonly("n_unitcells", &clexmonte::StateData::n_unitcells,
                    R"pbdoc(
          np.ndarray[np.int64] : The current state's supercell transformation \
          matrix.
          )pbdoc")
      .def_readonly("occ_location", &clexmonte::StateData::occ_location,
                    R"pbdoc(
          Optional[libcasm.monte.events.OccLocation] : The current state's occupant \
          location list. May be None.
          )pbdoc")
      .def_property_readonly(
          "convert",
          [](clexmonte::StateData &m) -> monte::Conversions const & {
            return *m.convert;
          },
          R"pbdoc(
          libcasm.monte.Conversions : Index conversions for the current state.
          )pbdoc")
      .def(
          "corr",
          [](clexmonte::StateData &m,
             std::string key) -> std::shared_ptr<clexulator::Correlations> {
            return m.corr.at(key);
          },
          R"pbdoc(
          Get a correlations calculator

          Parameters
          ----------
          key : str
              Basis set name

          Returns
          -------
          corr : libcasm.clexulator.Correlations
              The correlations calculator for `key`, set to calculate for
              `state`.
          )pbdoc",
          py::arg("key"))
      .def(
          "local_corr",
          [](clexmonte::StateData &m, std::string key)
              -> std::shared_ptr<clexulator::LocalCorrelations> {
            return m.local_corr.at(key);
          },
          R"pbdoc(
          Get a local correlations calculator

          Parameters
          ----------
          key : str
              Local basis set name

          Returns
          -------
          local_corr : libcasm.clexulator.LocalCorrelations
              The local correlations calculator for `key`, set to calculate for
              `state`.
          )pbdoc",
          py::arg("key"))
      .def(
          "clex",
          [](clexmonte::StateData &m,
             std::string key) -> std::shared_ptr<clexulator::ClusterExpansion> {
            return m.clex.at(key);
          },
          R"pbdoc(
          Get a cluster expansion calculator

          Parameters
          ----------
          key : str
              Cluster expansion name

          Returns
          -------
          clex : libcasm.clexulator.ClusterExpansion
              The cluster expansion calculator for `key`, set to calculate for
              `state`.
          )pbdoc",
          py::arg("key"))
      .def(
          "multiclex",
          [](clexmonte::StateData &m, std::string key)
              -> std::pair<std::shared_ptr<clexulator::MultiClusterExpansion>,
                           std::map<std::string, Index>> {
            return m.multiclex.at(key);
          },
          R"pbdoc(
          Get a multi-cluster expansion calculator

          Parameters
          ----------
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
          py::arg("key"))
      .def(
          "local_clex",
          [](clexmonte::StateData &m, std::string key)
              -> std::shared_ptr<clexulator::LocalClusterExpansion> {
            return m.local_clex.at(key);
          },
          R"pbdoc(
          Get a local cluster expansion

          Parameters
          ----------
          key : str
              Local cluster expansion name

          Returns
          -------
          local_clex : libcasm.clexulator.LocalClusterExpansion
              The local cluster expansion calculator for `key`, set to
              calculate for `state`.
          )pbdoc",
          py::arg("key"))
      .def(
          "local_multiclex",
          [](clexmonte::StateData &m, std::string key)
              -> std::pair<
                  std::shared_ptr<clexulator::MultiLocalClusterExpansion>,
                  std::map<std::string, Index>> {
            return m.local_multiclex.at(key);
          },
          R"pbdoc(
          Get a local multi-cluster expansion,

          Parameters
          ----------
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
          py::arg("key"))
      .def(
          "order_parameter",
          [](clexmonte::StateData &m,
             std::string key) -> std::shared_ptr<clexulator::OrderParameter> {
            return m.order_parameters.at(key);
          },
          R"pbdoc(
          Get an order parameter calculator

          Parameters
          ----------
          key : str
              The order parameter name

          Returns
          -------
          order_parameter : libcasm.clexulator.OrderParameter
              The order parameter calculator for `key`, set to calculate for
              the current state.
          )pbdoc",
          py::arg("key"))
      .def(
          "local_orbit_composition_calculator",
          [](clexmonte::StateData &m, std::string key)
              -> std::shared_ptr<clexmonte::LocalOrbitCompositionCalculator> {
            return m.local_orbit_composition_calculators.at(key);
          },
          R"pbdoc(
          Get a local orbit composition calculator

          Parameters
          ----------
          key : str
              The local orbit composition calculator name

          Returns
          -------
          local_orbit_composition_calculator : libcasm.clexmonte.LocalOrbitCompositionCalculator
              The local orbit composition calculator `key`, set to calculate for
              the current state.
          )pbdoc",
          py::arg("key"));

  py::class_<clexmonte::BaseMonteCalculator::kmc_data_type,
             std::shared_ptr<clexmonte::BaseMonteCalculator::kmc_data_type>>(
      m, "KineticsData",
      R"pbdoc(
      Access kinetics data used in a Monte Carlo method

      )pbdoc")
      .def(py::init<>(),
           R"pbdoc(
        .. rubric:: Constructor

        )pbdoc")
      .def_readonly("sampling_fixture_label",
                    &clexmonte::BaseMonteCalculator::kmc_data_type::
                        sampling_fixture_label,
                    R"pbdoc(
          str: The current sampling fixture label.

          This will be set to the current sampling fixture label at sampling time.
          )pbdoc")
      .def_property_readonly(
          "sampling_fixture",
          [](clexmonte::BaseMonteCalculator::kmc_data_type &self)
              -> sampling_fixture_type const & {
            if (self.sampling_fixture == nullptr) {
              throw std::runtime_error(
                  "Error in KineticsData.sampling_fixture: "
                  "This is not set until just before the first sample.");
            }
            return *self.sampling_fixture;
          },
          R"pbdoc(
          libcasm.clexmonte.run_management.SamplingFixture: A reference to the
          current sampling fixture.

          This will be set to the current sampling fixture at sampling time.
          )pbdoc")
      .def_readonly("total_rate",
                    &clexmonte::BaseMonteCalculator::kmc_data_type::total_rate,
                    R"pbdoc(
          float: This will be set to the total event rate at sampling time.
          )pbdoc")
      .def_readonly("time",
                    &clexmonte::BaseMonteCalculator::kmc_data_type::time,
                    R"pbdoc(
          float: Current simulation time when sampling occurs.

          For time-based sampling this will be equal to the sampling time
          and not determined by the time any event occurred.
          For count-based sampling, this will be equal to the time the n-th
          (by step or pass) event occurred, where n is the step or pass when
          sampling is due.
          )pbdoc")
      .def_readonly("prev_time",
                    &clexmonte::BaseMonteCalculator::kmc_data_type::prev_time,
                    R"pbdoc(
          dict[str, float]: Simulation time at last sample, by sampling fixture label.

          This will be set to store the time when the last sample
          was taken, with key equal to sampling fixture label. This is set to
          0.0 when the run begins.
          )pbdoc")
      .def_readonly(
          "unique_atom_id",
          &clexmonte::BaseMonteCalculator::kmc_data_type::unique_atom_id,
          R"pbdoc(
          list[int]: Unique atom ID for each atom currently in the system.

          The ID ``unique_atom_id[l]`` is the unique atom ID for the atom at the
          position given by ``atom_positions_cart[:,l]``.
          )pbdoc")
      .def_readonly(
          "prev_unique_atom_id",
          &clexmonte::BaseMonteCalculator::kmc_data_type::prev_unique_atom_id,
          R"pbdoc(
          dict[str, list[int]]: Unique atom ID for each atom at last sample, by
          sampling fixture label.

          The ID ``prev_unique_atom_id[label][l]`` is the unique atom ID for the
          atom at the position given by ``prev_atom_positions_cart[label][:,l]``.
          )pbdoc")
      .def_readonly(
          "atom_name_index_list",
          &clexmonte::BaseMonteCalculator::kmc_data_type::atom_name_index_list,
          R"pbdoc(
          list[int]: Set this to hold atom name indices for each column of the
          atom position matrices.

          When sampling, this will hold the atom name index for each column of
          the atom position matrices. The atom name index is an index into
          :func:`OccSystem.atom_name_list <libcasm.occ_events.OccSystem.atom_name_list>`.
          )pbdoc")
      .def_readonly(
          "atom_positions_cart",
          &clexmonte::BaseMonteCalculator::kmc_data_type::atom_positions_cart,
          R"pbdoc(
          np.ndarray[np.float[3,n_atoms]]: Current atom positions, as columns in
          Cartesian coordinates.

          Before a sample is taken, this will be updated to contain the current
          atom positions in Cartesian coordinates, with shape=(3, n_atoms).
          Sampling functions can use this and `prev_atom_positions_cart` to
          calculate displacements.
          )pbdoc")
      .def_readonly("prev_atom_positions_cart",
                    &clexmonte::BaseMonteCalculator::kmc_data_type::
                        prev_atom_positions_cart,
                    R"pbdoc(
          dict[str, np.ndarray[np.float[3,n_atoms]]]: Atom positions at last
          sample, as columns in Cartesian coordinates, by sampling fixture label.

          Before a sample is taken, this will be updated to contain the current
          atom positions in Cartesian coordinates, with shape=(3, n_atoms).
          Sampling functions can use this and `prev_atom_positions_cart` to
          calculate displacements.
          )pbdoc");

  py::class_<calculator_type, std::shared_ptr<calculator_type>>
      pyMonteCalculator(m, "MonteCalculator",
                        R"pbdoc(
      Interface for running Monte Carlo calculations
      )pbdoc");

  py::class_<potential_type>(m, "MontePotential",
                             R"pbdoc(
      Interface to potential calculators

      )pbdoc")
      .def(py::init<>(&make_potential),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        calculator : libcasm.clexmonte.MonteCalculator
            Monte Carlo calculator which implements the potential.
        state : libcasm.clexmonte.MonteCarloState
            The state to be calculated.
        )pbdoc",
           py::arg("calculator"), py::arg("state"))
      .def_property_readonly("state_data", &potential_type::state_data, R"pbdoc(
          libcasm.clexmonte.StateData: Data for the current
          state being calculated
          )pbdoc")
      .def("per_supercell", &potential_type::per_supercell, R"pbdoc(
          Calculate and return the potential per supercell, using current state data

          Returns
          -------
          value: float
              The potential per supercell for the current state
          )pbdoc")
      .def("per_unitcell", &potential_type::per_unitcell, R"pbdoc(
          Calculate and return the potential per unit cell, using current state data

          Returns
          -------
          value: float
              The potential per unit cell for the current state
          )pbdoc")
      .def("occ_delta_per_supercell", &potential_type::occ_delta_per_supercell,
           R"pbdoc(
          Calculate and return the change in potential per supercell, using current
          state data

          Parameters
          ----------
          linear_site_index: list[int]
              The linear site indices of the sites changing occupation
          new_occ: list[int]
              The new occupation indices on the sites.

          Returns
          -------
          value: float
              The change in potential per supercell from the current state
              to the state with the new occupation.
          )pbdoc",
           py::arg("linear_site_index"), py::arg("new_occ"));

  py::class_<clexmonte::PrimEventData>(m, "PrimEventData",
                                       R"pbdoc(
      Data common to all translationally equivalent events
      )pbdoc")
      .def_readonly("event_type_name",
                    &clexmonte::PrimEventData::event_type_name,
                    R"pbdoc(
          str: Event type name.
          )pbdoc")
      .def_readonly("equivalent_index",
                    &clexmonte::PrimEventData::equivalent_index,
                    R"pbdoc(
          int: Equivalent event index.
          )pbdoc")
      .def_readonly("is_forward", &clexmonte::PrimEventData::is_forward,
                    R"pbdoc(
          bool: Is forward trajectory (else reverse).
          )pbdoc")
      .def_readonly("prim_event_index",
                    &clexmonte::PrimEventData::prim_event_index,
                    R"pbdoc(
          int: Linear index for this prim event
          )pbdoc")
      .def_readonly("event", &clexmonte::PrimEventData::event,
                    R"pbdoc(
          libcasm.occ_events.OccEvent: Event definition
          )pbdoc")
      .def(
          "sites", [](clexmonte::PrimEventData &self) { return self.sites; },
          R"pbdoc(
          list[:class:`~libcasm.xtal.IntegralSiteCoordinate`]: Event sites,
          relative to origin unit cell
          )pbdoc")
      .def(
          "occ_init",
          [](clexmonte::PrimEventData &self) { return self.occ_init; },
          R"pbdoc(
          list[int]: Initial site occupation
          )pbdoc")
      .def(
          "occ_final",
          [](clexmonte::PrimEventData &self) { return self.occ_final; },
          R"pbdoc(
          list[int]: final site occupation
          )pbdoc")
      .def(
          "to_dict",
          [](clexmonte::PrimEventData &self,
             std::optional<std::reference_wrapper<occ_events::OccSystem const>>
                 system,
             bool include_cluster, bool include_cluster_occupation,
             bool include_event_invariants) -> nlohmann::json {
            jsonParser json;
            occ_events::OccEventOutputOptions opt;
            opt.include_cluster = include_cluster;
            opt.include_cluster_occupation = include_cluster_occupation;
            opt.include_event_invariants = include_event_invariants;
            to_json(self, json, system, opt);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the PrimEventData as a Python dict

          Parameters
          ----------
          event_system : Optional[libcasm.occ_events.OccSystem] = None
              A :class:`~libcasm.occ_events.OccSystem`. Providing `event_system`
              allows output of more event information, including occupant and
              atom names, cluster information, and symmetry information.

          include_cluster: bool = True
              If True, also include the cluster sites

          include_cluster_occupation: bool = True
              If True, also include the initial and final cluster occupation

          include_event_invariants: bool = True
              If True, also include event invariants: number of trajectories,
              number of each occupant type, and site distances

          Returns
          -------
          data : dict
              The PrimEventData as a Python dict
          )pbdoc",
          py::arg("event_system") = std::nullopt,
          py::arg("include_cluster") = true,
          py::arg("include_cluster_occupation") = true,
          py::arg("include_event_invariants") = true)
      .def("__repr__", [](clexmonte::PrimEventData &self) -> nlohmann::json {
        std::stringstream ss;
        jsonParser json;
        occ_events::OccEventOutputOptions opt;
        to_json(self, json, std::nullopt, opt);
        ss << json;
        return ss.str();
      });

  py::bind_vector<std::vector<clexmonte::PrimEventData>>(m, "PrimEventList",
                                                         R"pbdoc(
      PrimEventList is a list[:class:`PrimEventData`]-like object.
      )pbdoc");

  py::class_<clexmonte::EventID>(m, "EventID",
                                 R"pbdoc(
      Identifies an event via linear unit cell index in some supercell

      .. rubric:: Special Methods

      - Sort EventID using ``<``, ``<=``, ``>``, ``>=``, and compare
        using ``==`` and ``!=``
      - EventID may be copied with
        :func:`EventID.copy <libcasm.clexmonte.EventID.copy>`,
        `copy.copy`, or `copy.deepcopy`.
      - EventID is hashable and may be used as a key in a dict.

      )pbdoc")
      .def(py::init<Index, Index>(),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          prim_event_index: int
              Index specifying an event in
              :py:attr:`MonteEventData.prim_event_list <libcasm.clexmonte.MonteEventData.prim_event_list>`.
          unitcell_index: int:
              Linear unit cell index into a supercell, as determined by
              :class:`~libcasm.xtal.UnitCellIndexConverter`.
          )pbdoc",
           py::arg("prim_event_index"), py::arg("unitcell_index"))
      .def_readonly("prim_event_index", &clexmonte::EventID::prim_event_index,
                    R"pbdoc(
          int: Index specifying an event in
          :py:attr:`MonteEventData.prim_event_list <libcasm.clexmonte.MonteEventData.prim_event_list>`.
          )pbdoc")
      .def_readonly("unitcell_index", &clexmonte::EventID::unitcell_index,
                    R"pbdoc(
          int: Linear unit cell index into a supercell, as determined by
          :class:`~libcasm.xtal.UnitCellIndexConverter`.
          )pbdoc")
      .def(py::self < py::self, "Sorts EventID.")
      .def(py::self <= py::self, "Sorts EventID.")
      .def(py::self > py::self, "Sorts EventID.")
      .def(py::self >= py::self, "Sorts EventID.")
      .def(py::self == py::self, "Compare EventID.")
      .def(py::self != py::self, "Compare EventID.")
      .def("__hash__",
           [](clexmonte::EventID const &self) {
             return py::hash(
                 py::make_tuple(self.prim_event_index, self.unitcell_index));
           })
      .def(
          "copy",
          [](clexmonte::EventID const &self) {
            return clexmonte::EventID(self);
          },
          "Represent the EventID as a Python dict.")
      .def("__copy__",
           [](clexmonte::EventID const &self) {
             return clexmonte::EventID(self);
           })
      .def("__deepcopy__", [](clexmonte::EventID const &self,
                              py::dict) { return clexmonte::EventID(self); })
      .def(
          "to_dict",
          [](clexmonte::EventID const &self) {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          "Represent the EventID as a Python dict.")
      .def("__repr__",
           [](clexmonte::EventID const &self) {
             std::stringstream ss;
             jsonParser json;
             to_json(self, json);
             ss << json;
             return ss.str();
           })
      .def_static(
          "from_dict",
          [](nlohmann::json const &data) {
            // print errors and warnings to sys.stdout
            py::scoped_ostream_redirect redirect;
            jsonParser json{data};
            InputParser<clexmonte::EventID> event_id_parser(json);
            std::runtime_error error_if_invalid{
                "Error in libcasm.clexmonte.EventID.from_dict"};
            report_and_throw_if_invalid(event_id_parser, CASM::log(),
                                        error_if_invalid);
            return (*event_id_parser.value);
          },
          "Construct an EventID from a Python dict.");

  py::class_<clexmonte::EventData>(m, "EventData",
                                   R"pbdoc(
      Data particular to a single translationally distinct event

      Notes
      -----

      - EventData is obtained from
        :py:attr:`MonteEventData.event_data <libcasm.clexmonte.MonteEventData.event_data>`.
      - No constructor is provided

      )pbdoc")
      .def_readonly("event", &clexmonte::EventData::event,
                    R"pbdoc(
          libcasm.monte.events.OccEvent: Used to apply event and track occupants
          when the event is selected.
          )pbdoc")
      .def_readonly("unitcell_index", &clexmonte::EventData::unitcell_index,
                    R"pbdoc(
          int: Linear unit cell index into a supercell, as determined by
          :class:`~libcasm.xtal.UnitCellIndexConverter`.
          )pbdoc")
      .def(
          "to_dict",
          [](clexmonte::EventData const &self) {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          "Represent the EventData as a Python dict.")
      .def("__repr__", [](clexmonte::EventData const &self) {
        std::stringstream ss;
        jsonParser json;
        to_json(self, json);
        ss << json;
        return ss.str();
      });

  py::class_<clexmonte::EventState>(m, "EventState",
                                    R"pbdoc(
      Data calculated for a single event in a single state
      )pbdoc")
      .def_readonly("is_allowed", &clexmonte::EventState::is_allowed,
                    R"pbdoc(
          bool: True if event is allowed given current configuration; False otherwise.
          )pbdoc")
      .def_property_readonly(
          "formation_energy_delta_corr",
          [](clexmonte::EventState const &self) {
            if (self.formation_energy_delta_corr == nullptr) {
              throw std::runtime_error(
                  "Error in EventState.formation_energy_delta_corr: "
                  "not calculated.");
            }
            return *self.formation_energy_delta_corr;
          },
          R"pbdoc(
          numpy.ndarray[numpy.float[corr_size,]]: Change in formation energy
          correlations if event occurs.
          )pbdoc")
      .def_property_readonly(
          "local_corr",
          [](clexmonte::EventState const &self) {
            if (self.local_corr == nullptr) {
              throw std::runtime_error(
                  "Error in EventState.local_corr: "
                  "not calculated.");
            }
            return *self.local_corr;
          },
          R"pbdoc(
          numpy.ndarray[numpy.float[corr_size,]]: Local correlations for current
          event neighborhood.
          )pbdoc")
      .def_readonly("is_normal", &clexmonte::EventState::is_normal,
                    R"pbdoc(
          bool: An event is "normal" if `dE_activated` > 0.0 and
          `dE_activated` > `dE_final`.
          )pbdoc")
      .def_readonly("dE_final", &clexmonte::EventState::dE_final,
                    R"pbdoc(
          float: Final state energy, relative to initial state.
          )pbdoc")
      .def_readonly("Ekra", &clexmonte::EventState::Ekra,
                    R"pbdoc(
          float: KRA energy.
          )pbdoc")
      .def_readonly("dE_activated", &clexmonte::EventState::dE_activated,
                    R"pbdoc(
          float: Activated state energy, relative to initial state
          )pbdoc")
      .def_readonly("freq", &clexmonte::EventState::freq,
                    R"pbdoc(
          float: Attempt frequency
          )pbdoc")
      .def_readonly("rate", &clexmonte::EventState::rate,
                    R"pbdoc(
          float: Event rate
          )pbdoc")
      .def(
          "to_dict",
          [](clexmonte::EventState const &self) {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          "Represent the EventState as a Python dict.")
      .def("__repr__", [](clexmonte::EventState const &self) {
        std::stringstream ss;
        jsonParser json;
        to_json(self, json);
        ss << json;
        return ss.str();
      });

  py::class_<clexmonte::SelectedEvent,
             std::shared_ptr<clexmonte::SelectedEvent>>(m, "SelectedEvent",
                                                        R"pbdoc(
      Data structure holding the last selected event and its state
      )pbdoc")
      .def(py::init<>(),
           R"pbdoc(
          .. rubric:: Constructor

          Default constructor only.

          )pbdoc")
      .def_readonly("event_id", &clexmonte::SelectedEvent::event_id,
                    R"pbdoc(
          EventID: The event ID for the selected event.
          )pbdoc")
      .def_readonly("event_index", &clexmonte::SelectedEvent::event_index,
                    R"pbdoc(
          int: The index of the selected event in the event list (for
          `event_data_type` `"default"` or `"low_memory"`).
          )pbdoc")
      .def_readonly("total_rate", &clexmonte::SelectedEvent::total_rate,
                    R"pbdoc(
          float: The total rate when the event was selected.
          )pbdoc")
      .def_readonly("time_increment", &clexmonte::SelectedEvent::time_increment,
                    R"pbdoc(
          float: The time increment when the event occurred.
          )pbdoc")
      .def_readonly("prim_event_data",
                    &clexmonte::SelectedEvent::prim_event_data,
                    R"pbdoc(
          PrimEventData: Description of the selected event.
          )pbdoc")
      .def_readonly("event_data", &clexmonte::SelectedEvent::event_data,
                    R"pbdoc(
          EventData: Data for the selected event.
          )pbdoc")
      .def_readonly("event_state", &clexmonte::SelectedEvent::event_state,
                    R"pbdoc(
          Optional[EventState]: The calculated properties of the
          selected event.

          The kinetic Monte Carlo event lists may only contain the event rates
          and not store all the event properties. If
          `requires_event_state=True` for one or more selected event data
          functions, after the event is selected, the state of the selected
          event will be calculated and this will have a value; otherwise it
          will be None.
          )pbdoc");

  py::class_<clexmonte::MonteEventList>(m, "MonteEventList",
                                        R"pbdoc(
      Allows iterating over EventID

      .. rubric:: Special Methods

      - Use `len(event_list)` to get the number of events
      - Use ``for event_id in event_list:`` to iterate over the EventID

      )pbdoc")
      .def("__len__", &clexmonte::MonteEventList::size)
      .def("total_rate", &clexmonte::MonteEventList::total_rate)
      .def(
          "__iter__",
          [](clexmonte::MonteEventList const &self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<
              0, 1>() /* Essential: keep object alive while iterator exists */);

  py::class_<clexmonte::MonteEventData>(m, "MonteEventData",
                                        R"pbdoc(
      Interface to event data

      )pbdoc")
      .def(py::init<>(&make_event_data),
           R"pbdoc(
          .. rubric:: Constructor

          The constructor is equivalent to:

          .. code-block:: python

              if occ_location is None:
                  calculator.set_state_and_potential(state=state)
                  calculator.make_occ_location()
              else:
                  calculator.set_state_and_potential(
                      state=state,
                      occ_location=occ_location,
                  )
              calculator.set_event_data(engine=engine)
              event_data = calculator.event_data


          Parameters
          ----------
          calculator : libcasm.clexmonte.MonteCalculator
              Monte Carlo calculator which constructs events.
          state : libcasm.clexmonte.MonteCarloState
              The state events are constructed for. The calculator's state data
              and potential will be set to point to this state.
          engine: Optional[libcasm.monte.RandomNumberEngine] = None
              Optional random number engine to use. If None, one is constructed and
              seeded from std::random_device.
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state`. If None, an occupant location list owned by
              the calculator is constructed and initialized. The calculator's
              state data will be set to point to this occupant location list.

          )pbdoc",
           py::arg("calculator"), py::arg("state"), py::arg("engine") = nullptr,
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
      .def_property_readonly(
          "prim_event_list",
          [](clexmonte::MonteEventData &self)
              -> std::vector<clexmonte::PrimEventData> const & {
            return self.prim_event_list();
          },
          R"pbdoc(
          libcasm.clexmonte.PrimEventList: The translationally distinct
          instances of each event, including forward and reverse events separately,
          associated with origin primitive cell.
          )pbdoc")
      .def_property_readonly("event_list",
                             &clexmonte::MonteEventData::event_list,
                             R"pbdoc(
          MonteEventList: The current list of EventID.
          )pbdoc")
      .def("event_to_apply", &clexmonte::MonteEventData::event_to_apply,
           R"pbdoc(
          The event data structure that can used to apply the event to the current
          configuration's occupant location list.

          Parameters
          ----------
          id: libcasm.clexmonte.EventID
              The event ID for the event.

          Returns
          -------
          event: libcasm.monte.events.OccEvent
              The event data structure that can used to apply the event to the current
              configuration's occupant location list. The reference is valid until the
              next call to this function.
          )pbdoc")
      .def("event_rate", &clexmonte::MonteEventData::event_rate,
           R"pbdoc(
          Return the current rate for a specific event, as stored in the event list

          Parameters
          ----------
          id: libcasm.clexmonte.EventID
              The event ID for the event.

          Returns
          -------
          rate: float
              The current rate for the specified event, as stored in the event list.
          )pbdoc",
           py::arg("id"))
      .def("event_state", &clexmonte::MonteEventData::event_state,
           R"pbdoc(
          Calculate and return a reference to the EventState for a particular event
          in the current configuration

          Parameters
          ----------
          id: libcasm.clexmonte.EventID
              The event ID for the occuring event.

          Returns
          -------
          state: libcasm.clexmonte.EventState
              A reference to the EventState for a particular event in the current
              configuration. The reference is valid until the next call to
              this function.
          )pbdoc",
           py::arg("id"))
      .def("event_impact", &clexmonte::MonteEventData::event_impact,
           R"pbdoc(
          Return a list of EventID for the events that must be updated if a specified
          event occurs

          Parameters
          ----------
          id: libcasm.clexmonte.EventID
              The event ID for the occuring event.

          Returns
          -------
          impacted_events: list[libcasm.clexmonte.EventID]
              The EventID for events that must be updated.
          )pbdoc",
           py::arg("id"));

  py::class_<clexmonte::EventDataSummary,
             std::shared_ptr<clexmonte::EventDataSummary>>
      pyMonteEventDataSummary(m, "MonteEventDataSummary",
                              R"pbdoc(
      Summarizes MonteEventData

      )pbdoc");

  pyMonteCalculator
      .def(py::init<>(&make_monte_calculator),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          method : str
              Monte Carlo method name. The options are:

              - "semigrand_canonical": `Semi-grand canonical ensemble <todo>`_.
                Input states require `"temperature"` and `"param_chem_pot"`
                conditions.
              - "canonical": `Canonical ensemble <todo>`_.
                Input states require `"temperature"` and one of
                `"param_composition"` or `"mol_composition"` conditions.
              - TODO "lte": `Low-temperature expansion <todo>`_, for the
                semi-grand canonical ensemble
              - "kinetic": `Kinetic Monte Carlo <todo>`_. Input states require
                `"temperature"` and one of `"param_composition"` or
                `"mol_composition"` conditions.
              - TODO "flex": Allows a range of custom potentials, including
                composition and order parameter variance-constrained potentials,
                and correlation-matching potentials

          system : libcasm.clexmonte.System
              Cluster expansion model system data. The required data depends on
              the calculation method. See links under `method` for what system
              data is required for each method.

          params: Optional[dict] = None
              Monte Carlo calculation method parameters. Expected values
              depends on the calculation method. Options, with links to
              parameter documentation and examples, include:

              - "enumeration": `Save states <todo>`_ encountered during the
                calculation.

          )pbdoc",
           py::arg("method"), py::arg("system"),
           py::arg("params") = std::nullopt)
      .def(
          "make_default_sampling_fixture_params",
          [](std::shared_ptr<calculator_type> &self, std::string label,
             bool write_results, bool write_trajectory, bool write_observations,
             bool write_status, std::optional<std::string> output_dir,
             std::optional<std::string> log_file,
             double log_frequency_in_s) -> sampling_fixture_params_type {
            return self->make_default_sampling_fixture_params(
                self, label, write_results, write_trajectory,
                write_observations, write_status, output_dir, log_file,
                log_frequency_in_s);
          },
          R"pbdoc(
          Construct default sampling fixture parameters

          Notes
          -----

          By default:

          - Sampling occurs linearly, by pass, with period 1, for:

            - "clex.formation_energy": Formation energy, per unitcell
            - "mol_composition": Mol composition, :math:`\vec{n}`, per unitcell
            - "param_composition": Parametric composition, :math:\vec{x}`
            - "potential_energy": Potential energy, per unitcell
            - "order_parameter.<key>": Order parameter values (for all
              dof_spaces keys), and
            - "order_parameter.<key>.subspace_magnitudes": Magnitude of order
              parameter values in subspaces (for all dof_subspaces keys).

          - Analysis functions are evaluated for:

            - "heat_capacity",
            - "mol_susc" (excluding "canonical", "kinetic"),
            - "param_susc" (excluding "canonical", "kinetic"),
            - "mol_thermochem_susc" (excluding "canonical", "kinetic"), and
            - "param_thermochem_susc" (excluding "canonical", "kinetic").

          - Convergence of "potential_energy" is set to an
            absolute precision of 0.001, and "param_composition" to 0.001
            (excluding "canonical", "kinetic").
          - Completion is checked every 100 samples, starting with the 100-th.
          - No cutoffs are set.

          - Other standard sampling functions which are not included by default
            are:

            - "corr.<key>": Correlation values, for all basis functions (using
              basis sets key),
            - "clex.<key>": Cluster expansion value (using clex key),
            - "clex.<key>.sparse_corr": Correlations for cluster expansion
              basis functions with non-zero coefficients (using clex key),
            - "multiclex.<key>": Multi-cluster expansion value (using
              multiclex key),
            - "multiclex.<key>.sparse_corr" : Correlations for multi-cluster
              expansion basis functions with non-zero coefficients (using
              multiclex key),


          Parameters
          ----------
          label: str
              Label for the :class:`SamplingFixture`.
          write_results: bool = True
              If True, write results to summary file upon completion. If a
              results summary file already exists, the new results are appended.
          write_trajectory: bool = False
              If True, write the trajectory of Monte Carlo states when each
              sample taken to an output file. May be large.
          write_observations: bool = False
              If True, write a file with all individual sample observations.
              May be large.
          output_dir: Optional[str] = None
              Directory in which write results. If None, uses
              ``"output" / label``.
          write_status: bool = True
              If True, write log files with convergence status.
          log_file: str = Optional[str] = None
              Path to where a run status log file should be written with run
              information. If None, uses ``output_dir / "status.json"``.
          log_frequency_in_s: float = 600.0
              Minimum time between when the status log should be written, in
              seconds. The status log is only written after a sample is taken,
              so if the `sampling_params` are such that the time between
              samples is longer than `log_frequency_is_s` the status log will
              be written less frequently.

          Returns
          -------
          sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
              Default sampling fixture parameters for a semi-grand canonical
              Monte Carlo calculation.
          )pbdoc",
          py::arg("label"), py::arg("write_results") = true,
          py::arg("write_trajectory") = false,
          py::arg("write_observations") = false, py::arg("write_status") = true,
          py::arg("output_dir") = std::nullopt,
          py::arg("log_file") = std::nullopt,
          py::arg("log_frequency_in_s") = 600.0)
      .def(
          "make_sampling_fixture_params_from_dict",
          [](std::shared_ptr<calculator_type> &self, const nlohmann::json &data,
             std::string label) -> sampling_fixture_params_type {
            jsonParser json{data};
            bool time_sampling_allowed = false;

            InputParser<sampling_fixture_params_type> parser(
                json, label, self->sampling_functions,
                self->json_sampling_functions, self->analysis_functions,
                clexmonte::standard_results_io_methods(),
                time_sampling_allowed);
            std::runtime_error error_if_invalid{
                "Error in "
                "libcasm.clexmonte.MonteCalculator.sampling_fixture_params_"
                "from_dict"};
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            return std::move(*parser.value);
          },
          R"pbdoc(
          Construct sampling fixture parameters from Python dict

          Parameters
          ----------
          data: dict
              Python dict with sampling fixture parameters.
          label: str
              Label for the :class:`SamplingFixture`.

          Returns
          -------
          sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
              Sampling fixture parameters.
          )pbdoc",
          py::arg("data"), py::arg("label"))
      .def("set_state_and_potential", &calculator_type::set_state_and_potential,
           R"pbdoc(
          Set the current state and constructs the potential calculator

          Notes
          -----

          - Once called, it is possible to use the potential calculator,
            :py:attr:`~libcasm.clexmonte.MonteCalculator.potential`, outside of
            the `run` method.

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The input state.
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state` and it is used and updated during the run.
              The calculator's state data will be set to point to this occupant
              location list. If None, no occupant location list is stored. The
              occupant location list is not required for evaluating the potential.

          )pbdoc",
           py::arg("state"),
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
      .def("make_occ_location", &calculator_type::make_occ_location,
           R"pbdoc(
          Make and initialize an occupant location list for the current state

          Notes
          -----

          - Before calling, it is necessary to first set the state and potential,
            using :func:`~libcasm.clexmonte.MonteCalculator.set_state_and_potential`.
          - After calling, the calculator's state data will be updated with a
            pointer to the resulting occupant location list.

          Returns
          -------
          occ_location: libcasm.monte.events.OccLocation
              A current occupant location list initialized with the current
              state's occupation.

          )pbdoc")
      .def("set_event_data", &calculator_type::set_event_data,
           R"pbdoc(
          Set event data (includes calculating all rates), using current state data

          Notes
          -----

          - Before calling, it is necessary to first set the state and potential,
            using :func:`~libcasm.clexmonte.MonteCalculator.set_state_and_potential`,
            and make and set an occupant location list using
            :func:`~libcasm.clexmonte.MonteCalculator.make_occ_location`.


          Parameters
          ----------
          engine : Optional[:class:`~libcasm.monte.RandomNumberEngine`]
              A :class:`~libcasm.monte.RandomNumberEngine` to use for generating
              random numbers to select events and timesteps. If provided, the
              engine will be shared. If None, then a new
              :class:`~libcasm.monte.RandomNumberEngine` will be constructed and
              seeded using std::random_device.
          )pbdoc",
           py::arg("engine") = std::shared_ptr<engine_type>())
      .def("run", &monte_calculator_run,
           R"pbdoc(
          Perform a single run, evolving the input state

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The input state.
          run_manager: libcasm.clexmonte.RunManager
              Specifies sampling and convergence criteria and collects results
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state` and it is used and updated during the run.
              If None, a occupant location list is generated for the run.

          Returns
          -------
          run_manager: libcasm.clexmonte.RunManager
              The input `run_manager` with collected results.
          )pbdoc",
           py::arg("state"), py::arg("run_manager"),
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
      .def("run_fixture", &monte_calculator_run_fixture,
           R"pbdoc(
          Perform a single run, evolving the input state

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
            The input state.
          sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
              Specifies sampling and convergence criteria and collects results.
          engine: Optional[libcasm.monte.RandomNumberEngine] = None
              Optional random number engine to use. If None, one is constructed and
              seeded from std::random_device.
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state`. It is used and updated during the run.
              If None, an occupant location list is generated for the run.

          Returns
          -------
          sampling_fixture: libcasm.clexmonte.SamplingFixture
              A SamplingFixture with collected results.

          )pbdoc",
           py::arg("state"), py::arg("sampling_fixture_params"),
           py::arg("engine") = nullptr,
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
      .def_readwrite("sampling_functions", &calculator_type::sampling_functions,
                     R"pbdoc(
          libcasm.monte.StateSamplingFunctionMap: Sampling functions
          )pbdoc")
      .def_readwrite("json_sampling_functions",
                     &calculator_type::json_sampling_functions,
                     R"pbdoc(
          libcasm.monte.jsonStateSamplingFunctionMap: JSON sampling functions
          )pbdoc")
      .def_readwrite("analysis_functions", &calculator_type::analysis_functions,
                     R"pbdoc(
          libcasm.clexmonte.ResultsAnalysisFunctionMap: Results analysis functions
          )pbdoc")
      .def_readwrite("modifying_functions",
                     &calculator_type::modifying_functions,
                     R"pbdoc(
          libcasm.clexmonte.StateModifyingFunctionMap: State modifying functions
          )pbdoc")
      .def_property_readonly("selected_event", &calculator_type::selected_event,
                             R"pbdoc(
          Optional[libcasm.monte.sampling.SelectedEvent]: If applicable, will
          be set to provide information about the last selected event.

          If applicable for a particular calculator method, this will be set
          after each event selection to make the selected event information
          available for collection.

          If not applicable for a particular calculator method, this will be
          None.
          )pbdoc")
      .def_property_readonly("selected_event_functions",
                             &calculator_type::selected_event_functions,
                             R"pbdoc(
          Optional[libcasm.monte.sampling.SelectedEventFunctions]: Selected
          event data functions.

          If applicable for a particular calculator method, this will be
          constructed with standard functions and allows adding additional
          custom functions that can be called after each event is selected to
          collect selected event data, if requested by
          :func:`MonteCalculator.selected_event_function_params`.

          If not applicable for a particular calculator method, this will be
          None.
          )pbdoc")
      .def(
          "add_generic_function",
          [](calculator_type &self, std::string name, std::string description,
             bool requires_event_state, std::function<void()> function,
             std::function<bool()> has_value_function, Index order) {
            // -- Validation --
            if (function == nullptr) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_generic_function: "
                  "function=None");
            }

            std::shared_ptr<monte::SelectedEventFunctions> functions =
                self.selected_event_functions();
            if (!functions) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_generic_function: "
                  "selected_event_functions == None");
            }

            // -- End validation --

            if (has_value_function == nullptr) {
              has_value_function = []() { return true; };
            }

            monte::GenericSelectedEventFunction f(
                name, description, requires_event_state, function,
                has_value_function, order);

            functions->insert(f);
          },
          R"pbdoc(
          Add a function that will be evaluated after
          selecting an event and before collecting data about the event.

          This is a convenience method that:

          - constructs a
            :class:`~libcasm.monte.sampling.GenericSelectedEventFunction`
            with the given parameters, and
          - adds it to `self.selected_event_functions.generic_functions`.

          Parameters
          ----------
          name : str
              Name of the sampled quantity.
          description : str
              Description of the function.
          requires_event_state : bool
              If true, the function requires the event state of the selected
              event to be calculated.
          function : function
              A function with 0 arguments that returns a float. Typically this
              is a lambda function that has been given a reference or pointer to
              a Monte Carlo calculation object so that it can access the last
              selected event state and the Monte Carlo state before the event
              occurs.
          has_value_function : Optional[function] = None
              An optional function with 0 arguments that returns a bool
              indicating that the value of `function` should be collected.
              Default is to always return True. Typically this is a lambda
              function that has been given a reference or pointer to a Monte
              Carlo calculation object so that it can access the current
              selected event and determine if a value should be collected.
          order : int = 0
              Order in which the function is evaluated. Functions with lower
              order are evaluated first. If two functions have the same order,
              they are evaluated in lexicographical order by function name.
           )pbdoc",
          py::arg("name"), py::arg("description"),
          py::arg("requires_event_state"), py::arg("function"),
          py::arg("has_value_functions") = nullptr, py::arg("order") = 0)
      .def(
          "add_discrete_vector_int_function",
          [](calculator_type &self, std::string name, std::string description,
             std::vector<Index> shape, bool requires_event_state,
             std::function<Eigen::VectorXl()> function,
             std::function<bool()> has_value_function,
             std::optional<std::vector<std::string>> component_names,
             Index max_size,
             std::optional<std::vector<std::pair<Eigen::VectorXl, std::string>>>
                 value_labels) {
            // -- Validation --
            if (function == nullptr) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_discrete_vector_int_function: "
                  "function=None");
            }

            std::shared_ptr<monte::SelectedEventFunctions> functions =
                self.selected_event_functions();
            if (!functions) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_discrete_vector_int_function: "
                  "selected_event_functions == None");
            }

            // -- End validation --

            if (has_value_function == nullptr) {
              has_value_function = []() { return true; };
            }

            monte::DiscreteVectorIntHistogramFunction f(
                name, description, shape, component_names, requires_event_state,
                function, has_value_function, max_size);

            if (value_labels.has_value()) {
              for (auto const &pair : value_labels.value()) {
                f.value_labels->emplace(pair.first, pair.second);
              }
            }

            functions->insert(f);
          },
          R"pbdoc(
          Add a function for collecting integer-valued data after selecting an
          event.

          This is a convenience method that:

          - constructs a
            :class:`~libcasm.monte.sampling.VectorIntHistogramFunction`
            with the given parameters, and
          - adds it to `self.selected_event_functions.discrete_vector_int_functions`.

          Parameters
          ----------
          name : str
              Name of the sampled quantity.

          description : str
              Description of the function.

          shape : list[int]
              Shape of quantity, with column-major unrolling

              Scalar: [], Vector: [n], Matrix: [m, n], etc.

          requires_event_state : bool
              If true, the function requires the event state of the selected
              event to be calculated.

          function : function
              A function with 0 arguments that returns a integer-valued array.
              Typically this is a lambda function that has been given a
              reference or pointer to a Monte Carlo calculation object so that
              it can access the last selected event state and the Monte Carlo
              state before the event occurs.

          has_value_function : Optional[function] = None
              An optional function with 0 arguments that returns a bool
              indicating that the value of `function` should be collected.
              Default is to always return True. Typically this is a lambda
              function that has been given a reference or pointer to a Monte
              Carlo calculation object so that it can access the current
              selected event and determine if a value should be collected.

          component_names : Optional[list[str]] = None
              A name for each component of the resulting vector.

              Can be strings representing an indices (i.e "0", "1", "2", etc.) or can be a descriptive string (i.e. "Mg", "Va", "O", etc.). If None, indices for column-major ordering are used (i.e. "0,0", "1,0", ..., "m-1,n-1")

          max_size : int = 10000
              Maximum number of bins to create. If adding an additional data
              point would cause the number of bins to exceed `max_size`, the
              count / weight is instead added to the `out_of_range_count` of the
              :class:`~libcasm.monte.sampling.DiscreteVectorIntHistogram`.

          value_labels: Optional[list[tuple[np.ndarray, str]]]
              A list of tuples containing values and labels that if provided
              will be used to label the values in the output.

              For example, the standard function "selected_event.by_type"
              returns a value `[0]`, `[1]`, ... to indicate which type of event
              was selected and labels the output with the event type name
              using:

              .. code-block:: python

                  value_labels = [
                      (np.array([0]), "A_Va_1NN"),
                      (np.array([1]), "B_Va_1NN"),
                      ...
                  ]


          )pbdoc",
          py::arg("name"), py::arg("description"), py::arg("shape"),
          py::arg("requires_event_state"), py::arg("function"),
          py::arg("has_value_function") = nullptr,
          py::arg("component_names") = std::nullopt,
          py::arg("max_size") = 10000, py::arg("value_labels") = std::nullopt)
      .def(
          "add_discrete_vector_float_function",
          [](calculator_type &self, std::string name, std::string description,
             std::vector<Index> shape, bool requires_event_state,
             std::function<Eigen::VectorXd()> function,
             std::function<bool()> has_value_function,
             std::optional<std::vector<std::string>> component_names,
             Index max_size, double tol,
             std::optional<std::vector<std::pair<Eigen::VectorXd, std::string>>>
                 value_labels) {
            // -- Validation --
            if (function == nullptr) {
              throw std::runtime_error(
                  "Error in "
                  "MonteCalculator.add_discrete_vector_float_function: "
                  "function=None");
            }

            std::shared_ptr<monte::SelectedEventFunctions> functions =
                self.selected_event_functions();
            if (!functions) {
              throw std::runtime_error(
                  "Error in "
                  "MonteCalculator.add_discrete_vector_float_function: "
                  "selected_event_functions == None");
            }

            // -- End validation --

            if (has_value_function == nullptr) {
              has_value_function = []() { return true; };
            }

            monte::DiscreteVectorFloatHistogramFunction f(
                name, description, shape, component_names, requires_event_state,
                function, has_value_function, max_size, tol);

            if (value_labels.has_value()) {
              for (auto const &pair : value_labels.value()) {
                f.value_labels->emplace(pair.first, pair.second);
              }
            }

            functions->insert(f);
          },
          R"pbdoc(
          Add a function for collecting discrete floating-valued data after
          selecting an event.

          This is a convenience method that:

          - constructs a
            :class:`~libcasm.monte.sampling.VectorFloatHistogramFunction`
            with the given parameters, and
          - adds it to `self.selected_event_functions.discrete_vector_float_functions`.

          Parameters
          ----------
          name : str
              Name of the sampled quantity.

          description : str
              Description of the function.

          shape : list[int]
              Shape of quantity, with column-major unrolling

              Scalar: [], Vector: [n], Matrix: [m, n], etc.

          requires_event_state : bool
              If true, the function requires the event state of the selected
              event to be calculated.

          function : function
              A function with 0 arguments that returns a float-valued array.
              Typically this is a lambda function that has been given a
              reference or pointer to a Monte Carlo calculation object so that
              it can access the last selected event state and the Monte Carlo
              state before the event occurs.

          has_value_function : Optional[function] = None
              An optional function with 0 arguments that returns a bool
              indicating that the value of `function` should be collected.
              Default is to always return True. Typically this is a lambda
              function that has been given a reference or pointer to a Monte
              Carlo calculation object so that it can access the current
              selected event and determine if a value should be collected.

          component_names : Optional[list[str]] = None
              A name for each component of the resulting vector.

              Can be strings representing an indices (i.e "0", "1", "2", etc.) or can be a descriptive string (i.e. "Mg", "Va", "O", etc.). If None, indices for column-major ordering are used (i.e. "0,0", "1,0", ..., "m-1,n-1")

          max_size : int = 10000
              Maximum number of bins to create. If adding an additional data
              point would cause the number of bins to exceed `max_size`, the
              count / weight is instead added to the `out_of_range_count` of the
              :class:`~libcasm.monte.sampling.DiscreteVectorIntHistogram`.

          tol : float = :data:`~libcasm.casmglobal.TOL`
              Tolerance for floating point comparisons used when determining
              counts for the histogram.

          value_labels: Optional[list[tuple[np.ndarray, str]]]
              A list of tuples containing values and labels that if provided
              will be used to label the values in the output.

              For example:

              .. code-block:: python

                  value_labels = [
                      (np.array([-1.0]), "negative"),
                      (np.array([1.0]), "positive"),
                      ...
                  ]


          )pbdoc",
          py::arg("name"), py::arg("description"), py::arg("shape"),
          py::arg("requires_event_state"), py::arg("function"),
          py::arg("has_value_function") = nullptr,
          py::arg("component_names") = std::nullopt,
          py::arg("max_size") = 10000, py::arg("tol") = CASM::TOL,
          py::arg("value_labels") = std::nullopt)
      .def(
          "add_partitioned_histogram_function",
          [](calculator_type &self, std::string name, std::string description,
             bool requires_event_state, std::function<double()> function,
             std::string partition_type, bool is_log, double initial_begin,
             double bin_width, Index max_size) {
            // -- Validation --
            if (function == nullptr) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_histogram_function: "
                  "function=None");
            }

            std::shared_ptr<monte::SelectedEventFunctions> functions =
                self.selected_event_functions();
            if (!functions) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_histogram_function: "
                  "selected_event_functions == None");
            }

            std::shared_ptr<clexmonte::SelectedEvent> selected_event =
                self.selected_event();
            if (!selected_event) {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_histogram_function: "
                  "selected_event == None");
            }

            // -- End validation --

            // This makes:
            // - a lookup table for prim_event_index -> partition index
            // - a vector of partition names
            clexmonte::SelectedEventInfo info(
                self.event_data().prim_event_list());
            if (partition_type == "by_type") {
              info.make_indices_by_type();
            } else if (partition_type == "by_equivalent_index") {
              info.make_indices_by_equivalent_index();
            } else if (partition_type == "by_equivalent_index_and_direction") {
              info.make_indices_by_equivalent_index_and_direction();
            } else {
              throw std::runtime_error(
                  "Error in MonteCalculator.add_histogram_function: "
                  "invalid partition_type=" +
                  partition_type);
            }

            // This is the lookup table, as a shared_ptr which can be
            // captured by the lambda function `get_partition`
            std::shared_ptr<std::vector<Index>> prim_event_index_to_index =
                info.prim_event_index_to_index;
            // The `selected_event` is a shared_ptr<SelectedEvent> which
            // will be updated by the MonteCalculator after each event is
            // selected.
            auto get_partition = [prim_event_index_to_index, selected_event]() {
              return prim_event_index_to_index->at(
                  selected_event->prim_event_data->prim_event_index);
            };

            monte::PartitionedHistogramFunction<double> f(
                name, description, requires_event_state, function,
                info.partition_names, get_partition, is_log, initial_begin,
                bin_width, max_size);

            functions->insert(f);
          },
          R"pbdoc(
          Add a custom function for collecting continuous scalar quantities of
          the selected events in 1d histograms, partitioned by (i) event type,
          or (ii) event type and equivalent index, or (iii) event type,
          equivalent index, and hop direction.

          This is a convenience method that:

          - constructs a lookup table (`prim_event_index` -> `partition_index`)
            to partition collected data by the type of selected event,
          - generates the corresponding list of partition names,
          - constructs a
            :class:`~libcasm.monte.sampling.PartitionedHistogramFunction`
            with the partitions and given parameters, and
          - adds it to `self.selected_event_functions.continuous_1d_functions`.

          Parameters
          ----------
          name : str
              Name of the sampled quantity.
          description : str
              Description of the function.
          requires_event_state : bool
              If true, the function requires the event state of the selected
              event to be calculated.
          function : function
              A function with 0 arguments that returns a float. Typically this
              is a lambda function that has been given a reference or pointer to
              a Monte Carlo calculation object so that it can access the last
              selected event state and the Monte Carlo state before the event
              occurs.
          partition_type : str
              Determines the type of partitioning to use for the histogram.
              Options are:

              - "by_type": Partition by event type
              - "by_equivalent_index": Partition by event type and equivalent
                index.
              - "by_equivalent_index_and_direction": Partition by event
                type, equivalent index, and direction

          is_log : bool = False
              True if bin coordinate spacing is log-scaled; False otherwise.
          initial_begin : float = 0.0
              Initial `begin` coordinate, specifying the beginning of the range
              for the first bin. The bin number for a particular value is
              calculated as `(value - begin) / bin_width`, so the range for
              bin `i` is [begin, begin + i*bin_width). Coordinates are adjusted
              to fit the data encountered by starting `begin` at
              `initial_begin` and adjusting it as necessary by multiples of
              `bin_width`.
          bin_width : float = 1.0
              Bin width.
          max_size : int = 10000
              Maximum number of bins to create. If adding an additional data
              point would cause the number of bins to exceed `max_size`, the
              count / weight is instead added to the `out_of_range_count` of the
              :class:`~libcasm.monte.sampling.PartitionedHistogram1D`.

          )pbdoc",
          py::arg("name"), py::arg("description"),
          py::arg("requires_event_state"), py::arg("function"),
          py::arg("partition_type"), py::arg("is_log") = false,
          py::arg("initial_begin") = 0.0, py::arg("bin_width") = 1.0,
          py::arg("max_size") = 10000)
      .def_property(
          "selected_event_function_params",
          &calculator_type::selected_event_function_params,
          [](calculator_type &self,
             std::shared_ptr<monte::SelectedEventFunctionParams>
                 selected_event_function_params) {
            self.set_selected_event_function_params(
                selected_event_function_params);
          },
          R"pbdoc(
          Optional[libcasm.monte.sampling.SelectedEventFunctionParams]: Selected event
          data collection parameters.

          If applicable for a particular calculator method, this may be set
          directly or read from the MonteCalculator constructor `params` input,
          if that is implemented by the calculator method. If it exists when a
          run begins, selected event data will be collected accordingly.

          If not applicable for a particular calculator method, this will be
          None.

          )pbdoc")
      .def(
          "collect_hop_correlations",
          [](calculator_type &self, Index jumps_per_position_sample,
             Index max_n_position_samples, bool output_incomplete_samples,
             bool stop_run_when_complete) {
            auto params_ptr = self.selected_event_function_params();
            if (params_ptr == nullptr) {
              self.set_selected_event_function_params(
                  std::make_shared<monte::SelectedEventFunctionParams>());
              params_ptr = self.selected_event_function_params();
            }
            params_ptr->correlations_data_params =
                monte::CorrelationsDataParams(
                    {jumps_per_position_sample, max_n_position_samples,
                     output_incomplete_samples, stop_run_when_complete});
            return self;
          },
          R"pbdoc(
          Update :py:attr:`selected_event_function_params` to collect hop
          correlations data

          Parameters
          ----------
          jumps_per_position_sample : int = 1
              Every `jumps_per_position_sample` steps of an individual atom,
              its position will be stored in Cartesian coordinates (as if
              periodic boundaries did not exist).
          max_n_position_samples : int = 100
              The maximum number of positions to store for each atom.
          output_incomplete_samples : bool = False
              If false, when representing this object as a Python dict, only
              output data for the number of samples for which all atoms have
              jumped the necessary number of times. If true, output matrices
              with 0.0 values for atoms that have not jumped enough times to be
              sampled.
          stop_run_when_complete : bool = False
              If true, stop the run when all atoms have jumped the necessary number
              of times to be sampled.

          Returns
          -------
          self: libcasm.clexmonte.SelectedEventFunctionParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("jumps_per_position_sample") = 1,
          py::arg("max_n_position_samples") = 100,
          py::arg("output_incomplete_samples") = false,
          py::arg("stop_run_when_complete") = false)
      .def(
          "do_not_collect_hop_correlations",
          [](calculator_type &self) {
            auto params_ptr = self.selected_event_function_params();
            if (params_ptr) {
              params_ptr->correlations_data_params = std::nullopt;
            }
            return self;
          },
          R"pbdoc(
          Update :py:attr:`selected_event_function_params` to not collect hop correlations data

          Returns
          -------
          self: libcasm.clexmonte.SelectedEventFunctionParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc")
      .def(
          "evaluate",
          [](calculator_type &self, std::string name,
             std::optional<Index> order = std::nullopt) {
            auto params_ptr = self.selected_event_function_params();
            if (params_ptr == nullptr) {
              self.set_selected_event_function_params(
                  std::make_shared<monte::SelectedEventFunctionParams>());
              params_ptr = self.selected_event_function_params();
            }
            params_ptr->evaluate(name, order);
            return self;
          },
          R"pbdoc(
          Update :py:attr:`selected_event_function_params` to add the name of a
          generic function that will be evaluated for each selected
          event, along with optional order of evaluation.

          Parameters
          ----------
          name : str
              The name of a Selected event function to be added to
              `self.selected_event_function_params.function_names`. These should
              be keys in one of the dictionaries in
              :func:`MonteCalculator.selected_event_functions`.
          order : Optional[int] = None
              The order in which generic selected event functions are evaluated.
              Functions with lower order are evaluated first. If two functions
              have the same order, they are evaluated in lexicographical order
              by function name.

          Returns
          -------
          self: libcasm.clexmonte.SelectedEventFunctionParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"), py::arg("order") = std::nullopt)
      .def(
          "collect",
          [](calculator_type &self, std::string name,
             std::optional<double> tol = std::nullopt,
             std::optional<double> bin_width = std::nullopt,
             std::optional<double> initial_begin = std::nullopt,
             std::optional<std::string> spacing = std::nullopt,
             std::optional<int> max_size = std::nullopt) {
            auto params_ptr = self.selected_event_function_params();
            if (params_ptr == nullptr) {
              self.set_selected_event_function_params(
                  std::make_shared<monte::SelectedEventFunctionParams>());
              params_ptr = self.selected_event_function_params();
            }
            params_ptr->collect(name, tol, bin_width, initial_begin, spacing,
                                max_size);
            return self;
          },
          R"pbdoc(
          Update :py:attr:`selected_event_function_params` to add the name of a
          function that will be evaluated to collect data for each selected
          event, along with optional custom settings.

          Parameters
          ----------
          name : str
              The name of a Selected event function to be added to
              `self.selected_event_function_params.function_names`. These should
              be keys in one of the dictionaries in
              :func:`MonteCalculator.selected_event_functions`.
          tol : Optional[float] = None
              The tolerance for comparing values, applicable to
              discrete floating point valued functions. If None, the
              function's default tolerance value is used.
          bin_width : Optional[float] = None
              The tolerance for comparing values, applicable to
              continuous valued functions. If None, the function's default bin
              width is used.
          initial_begin : Optional[float] = None
              The initial value for the first bin, applicable to continuous
              valued functions. If None, the function's default initial begin
              value is used.
          spacing : Optional[str] = None
              The spacing of the bins, applicable to continuous valued
              functions. If None, the function's default spacing is used.
              Options are "log" or "linear".
          max_size : Optional[int] = None
              The maximum number of bins to store, applicable to all functions.
              If None, the function's default maximum number of bins is used.

          Returns
          -------
          self: libcasm.clexmonte.SelectedEventFunctionParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"), py::arg("tol") = std::nullopt,
          py::arg("bin_width") = std::nullopt,
          py::arg("initial_begin") = std::nullopt,
          py::arg("spacing") = std::nullopt, py::arg("max_size") = std::nullopt)
      .def(
          "do_not_collect",
          [](calculator_type &self, std::string name) {
            auto params_ptr = self.selected_event_function_params();
            if (params_ptr) {
              params_ptr->do_not_collect(name);
            }
            return self;
          },
          R"pbdoc(
          Update :py:attr:`selected_event_function_params` to remove the name of a
          function that will be evaluated to collect data for each selected
          event, and to remove all custom settings.
          )pbdoc")
      .def_property_readonly("selected_event_data",
                             &calculator_type::selected_event_data,
                             R"pbdoc(
          Optional[libcasm.monte.sampling.SelectedEventData]: Selected event
          data

          If applicable for a particular calculator method, and requested by
          setting :func:`MonteCalculator.selected_event_function_params`, this will
          store selected event data.

          If not applicable for a particular calculator method, this will be
          None.
          )pbdoc")
      .def_property_readonly("name", &calculator_type::calculator_name,
                             R"pbdoc(
          str : Calculator name.
          )pbdoc")
      .def_property_readonly("system", &calculator_type::system, R"pbdoc(
          System : System data.
          )pbdoc")
      .def_property_readonly("time_sampling_allowed",
                             &calculator_type::time_sampling_allowed, R"pbdoc(
          bool : True if this calculation allows time-based sampling; \
          False otherwise.
          )pbdoc")
      .def_property_readonly("state_data", &calculator_type::state_data,
                             R"pbdoc(
          StateData : The current state data.
          )pbdoc")
      .def_property_readonly("potential", &calculator_type::potential, R"pbdoc(
          MontePotential : The potential calculator for the current state.
          )pbdoc")
      .def_property_readonly("run_manager", &calculator_type::run_manager,
                             R"pbdoc(
          RunManager : The current :class:`libcasm.clexmonte.RunManager`.
          )pbdoc")
      .def_property_readonly("event_data", &calculator_type::event_data,
                             R"pbdoc(
          MonteEventData : The current event data.
          )pbdoc")
      .def_property_readonly("kinetics_data", &calculator_type::kmc_data,
                             R"pbdoc(
          KineticsData : The current kinetics data.
          )pbdoc")
      .def(
          "event_data_summary",
          [](calculator_type &calculator, double energy_bin_width,
             double freq_bin_width, double rate_bin_width) {
            return std::make_shared<clexmonte::EventDataSummary>(
                calculator.state_data(), calculator.event_data(),
                energy_bin_width, freq_bin_width, rate_bin_width);
          },
          R"pbdoc(
          Construct :class:`~libcasm.clexmonte.MonteEventDataSummary` for the
          current event data

          Parameters
          ----------
          energy_bin_width : float = 0.1
              The bin width for energy histograms (eV).

          freq_bin_width : float = 0.1
              The bin width for frequency histograms (log10(1/s)).

          rate_bin_width : float = 0.1
              The bin width for rate histograms (log10(1/s)).

          Returns
          -------
          event_data_summary: libcasm.clexmonte.MonteEventDataSummary
              The event data summary. Includes:

              - "n_unitcells": The number of unit cells in the calculation
                supercell.
              - "n_events": The number of events in the current
                state in total, by event type, and by event type and
                equivalent event index.
              - "n_events_with_no_barrier": The number of events in the
                current state with no barrier (using the current model
                parameters) in total, by event type, and by event type and
                equivalent event index.
              - "event_list_size": The total event list size. This may be
                larger than the number of events in the current state if the
                event list method saves slots for events that are possible
                but not allowed in the current state (these are given rate 0.0
                and not chosen).
              - "rate": The sum of the rate of all events in total, events by
                event type, and events by event type and equivalent event index.
              - "mean_time_increment": The mean time increment until the next
                event occurs in total (1/total_rate), by event type, and by
                event type and equivalent event index.
              - "memory_used": The total memory used by the current process, as
                a str with units.
              - "memory_used_MiB": The total memory used by the current process
                as a float in MiB.
              - "impact_neighborhood": The number of sites where a change in DoF
                value triggers an update of the event rate, by event type.
                The number of sites that trigger an update due to the formation
                energy cluster expansion, the kra local cluster expansion, and
                the attempt frequency local cluster expansion are given
                individually, and the total number of sites, which is the union
                of the three, is also given. For example, if
                "impact_neighborhood/A_Va_1NN/total" has the value 20, then
                there are 20 sites where a change in the DoF values triggers an
                update of the `A_Va_1NN` event rate.
              - "impact_table": The number of events for which an
                update is triggered if a specified event occurs, by occurring
                event type (or event type and equivalent index), and by
                impacted event type (or event type and equivalent index)
                averaged over all events in the current event list (which may
                include both allowed and not allowed events, depending on the
                event list method), as tables, where the rows correspond to the
                occuring event type and the columns correspond to the impacted
                event type. The `"type"` and `"equiv"` attributes give the
                order of event type (or event type and equivalent index), for
                the rows and columns of the tables.

          )pbdoc",
          py::arg("energy_bin_width") = 0.1, py::arg("freq_bin_width") = 0.1,
          py::arg("rate_bin_width") = 0.1);

  m.def("make_custom_monte_calculator", &make_custom_monte_calculator, R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          system : libcasm.clexmonte.System
              Cluster expansion model system data. The required data depends on
              the calculation method.

          source: str
              Path to a MonteCalculator source file implementing a custom Monte
              Carlo method to use instead of a standard implementation.

          params: Optional[dict] = None
              Monte Carlo calculation method parameters. Expected values
              depends on the calculation method. Options, with links to
              parameter documentation and examples, include:

              - "enumeration": `Save states <todo>`_ encountered during the
                calculation.

          compile_options: Optional[str] = None
              Options used to compile the MonteCalculator source file, if it is not yet
              compiled. Example: "g++ -O3 -Wall -fPIC --std=c++17 -I/path/to/include".
              The default values can be configured with:

                  CASM_CXX:
                      Set compiler; default="g++"
                  CASM_CXXFLAGS:
                      Set compiler flags; default="-O3 -Wall -fPIC --std=c++17"
                  CASM_INCLUDEDIR:
                      Set include search path, overriding CASM_PREFIX
                  CASM_PREFIX:
                      Set include search path to -I$CASM_PREFIX/include; default
                      tries to find "ccasm" or "casm" executables on PATH and
                      checks relative locations

          so_options: Optional[str] = None
              Options used to compile the MonteCalculator shared object file, if it is not
              yet compiled. Example: "g++ -shared -L/path/to/lib -lcasm_clexmonte "

              The default values can be configured with:

                  CASM_CXX:
                      Set compiler; default="g++"
                  CASM_SOFLAGS:
                      Set shared object compilation flags; default="-shared"
                  CASM_LIBDIR:
                      Set link search path, overriding CASM_PREFIX
                  CASM_PREFIX:
                      Set include search path to -L$CASM_PREFIX/lib; default
                      tries to find "ccasm" or "casm" executables on PATH and
                      checks relative locations

          search_path: Optional[list[str]] = None
              An optional search path for the `source` file.

          )pbdoc",
        py::arg("system"), py::arg("source"), py::arg("params") = std::nullopt,
        py::arg("compile_options") = std::nullopt,
        py::arg("so_options") = std::nullopt,
        py::arg("search_path") = std::nullopt);

  pyMonteEventDataSummary
      .def(py::init<>(&make_event_data_summary),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          calculator : libcasm.clexmonte.MonteCalculator
              Monte Carlo calculator with event data already generated.
          energy_bin_width : float = 0.1
              The bin width for energy histograms (eV).
          freq_bin_width : float = 0.1
              The bin width for frequency histograms (log10(1/s)).
          rate_bin_width : float = 0.1
              The bin width for rate histograms (log10(1/s)).
          )pbdoc",
           py::arg("calculator"), py::arg("energy_bin_width") = 0.1,
           py::arg("freq_bin_width") = 0.1, py::arg("rate_bin_width") = 0.1)
      .def("__repr__",
           [](clexmonte::EventDataSummary const &event_data_summary) {
             std::stringstream ss;
             jsonParser json;
             auto const &x = event_data_summary;
             json["rate_total"]["total"] = x.total_rate;
             json["n_events_total"]["total"] = x.n_events_allowed;
             json["memory_used"] = convert_size(x.resident_bytes_used);
             ss << json;
             return ss.str();
           })
      .def("__str__",
           [](clexmonte::EventDataSummary const &event_data_summary) {
             std::stringstream ss;
             print(ss, event_data_summary);
             return ss.str();
           })
      .def(
          "to_dict",
          [](clexmonte::EventDataSummary const &event_data_summary)
              -> nlohmann::json {
            jsonParser json;
            to_json(event_data_summary, json);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Return MonteEventDataSummary as a dict

          Returns
          -------
          data: dict
              The event data summary as a dict. Includes:

              - "n_unitcells": The number of unit cells in the calculation
                supercell.
              - "n_events": The number of events in the current
                state in total, by event type, and by event type and
                equivalent event index.
              - "n_events_with_no_barrier": The number of events in the
                current state with no barrier (using the current model
                parameters) in total, by event type, and by event type and
                equivalent event index.
              - "event_list_size": The total event list size. This may be
                larger than the number of events in the current state if the
                event list method saves slots for events that are possible
                but not allowed in the current state (these are given rate 0.0
                and not chosen).
              - "rate": The sum of the rate of all events in total, events by
                event type, and events by event type and equivalent event index.
              - "mean_time_increment": The mean time increment until the next
                event occurs in total (1/total_rate), by event type, and by
                event type and equivalent event index.
              - "memory_used": The total memory used by the current process, as
                a str with units.
              - "memory_used_MiB": The total memory used by the current process
                as a float in MiB.
              - "impact_neighborhood": The number of sites where a change in DoF
                value triggers an update of the event rate, by event type.
                The number of sites that trigger an update due to the formation
                energy cluster expansion, the kra local cluster expansion, and
                the attempt frequency local cluster expansion are given
                individually, and the total number of sites, which is the union
                of the three, is also given. For example, if
                "impact_neighborhood/A_Va_1NN/total" has the value 20, then
                there are 20 sites where a change in the DoF values triggers an
                update of the `A_Va_1NN` event rate.
              - "impact_table": The number of events for which an
                update is triggered if a specified event occurs, by occurring
                event type (or event type and equivalent index), and by
                impacted event type (or event type and equivalent index)
                averaged over all allowed events, as tables, where the rows
                correspond to the occurring event type and the columns
                correspond to the impacted event type. The `"type"` and
                `"equiv"` attributes give the order of event type (or event
                type and equivalent index), for the rows and columns of the
                tables.

          )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
