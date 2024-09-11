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

// clexmonte/semigrand_canonical
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/clexmonte/monte_calculator/io/json/MonteCalculator_json_io.hh"
#include "casm/clexmonte/run/StateModifyingFunction.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
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
  monte::OccLocation *occ_location = nullptr;
  calculator->set_state_and_potential(state, nullptr);
  return calculator->potential();
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
  // Need to check for an OccLocation
  std::unique_ptr<monte::OccLocation> tmp;
  make_temporary_if_necessary(state, occ_location, tmp, self);

  // run
  self.run(state, *occ_location, *run_manager);
  return run_manager;
}

std::shared_ptr<sampling_fixture_type> monte_calculator_run_fixture(
    calculator_type &self, state_type &state,
    sampling_fixture_params_type &sampling_fixture_params,
    std::shared_ptr<engine_type> engine, monte::OccLocation *occ_location) {
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
              -> std::shared_ptr<clexulator::MultiClusterExpansion> {
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
              -> std::shared_ptr<clexulator::MultiLocalClusterExpansion> {
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
              If None, no occupant location list is stored. The occupant
              location list is not required for evaluating the potential.

          )pbdoc",
           py::arg("state"),
           py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr))
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
      .def_property_readonly("kinetics_data", &calculator_type::kmc_data,
                             R"pbdoc(
          KineticsData : The current kinetics data.
          )pbdoc");

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

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
