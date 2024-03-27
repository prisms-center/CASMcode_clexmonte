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
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/semigrand_canonical/calculator_impl.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexmonte/semigrand_canonical/json_io.hh"
#include "casm/clexmonte/semigrand_canonical/potential.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/monte/sampling/RequestedPrecisionConstructor.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;
namespace sgc = clexmonte::semigrand_canonical;

// used for libcasm.clexmonte:
typedef std::mt19937_64 engine_type;
typedef monte::RandomNumberGenerator<engine_type> generator_type;
typedef sgc::SemiGrandCanonical<engine_type> calculator_type;
typedef sgc::SemiGrandCanonicalPotential potential_type;
typedef sgc::SemiGrandCanonicalConditions conditions_type;
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

std::shared_ptr<calculator_type> make_calculator(
    std::shared_ptr<system_type> system) {
  std::shared_ptr<calculator_type> calculator =
      std::make_shared<calculator_type>(system);
  calculator->sampling_functions =
      calculator_type::standard_sampling_functions(calculator);
  calculator->json_sampling_functions =
      calculator_type::standard_json_sampling_functions(calculator);
  calculator->analysis_functions =
      calculator_type::standard_analysis_functions(calculator);
  return calculator;
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MAKE_OPAQUE(CASM::monte::SamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonSamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::StateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonStateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASMpy::analysis_function_map_type);

PYBIND11_MODULE(_clexmonte_semigrand_canonical, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
    Cluster expansion semi-grand canonical Monte Carlo
    )pbdoc";
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.events");
  py::module::import("libcasm.monte.sampling");
  py::module::import("libcasm.clexmonte");

  // TODO: Document methods
  py::class_<conditions_type, std::shared_ptr<conditions_type>>(
      m, "SemiGrandCanonicalConditions",
      R"pbdoc(
      Semi-grand canonical Monte Carlo thermodynamic conditions
      )pbdoc")
      .def(py::init<composition::CompositionConverter const &, double>(),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        composition_converter : libcasm.composition.CompositionConverter
            Cluster expansion model system data.
        temperature_is_zero_tol: float = 1e-10
        )pbdoc",
           py::arg("composition_converter"),
           py::arg("temperature_is_zero_tol") = 1e-10)
      .def_readonly("temperature_is_zero_tol",
                    &conditions_type::temperature_is_zero_tol, R"pbdoc(
          float: Tolerance used to check if :math`T=0`, in which case `beta` \
          is set to infinity.
          )pbdoc")
      .def_readonly("temperature", &conditions_type::temperature, R"pbdoc(
          float: The temperature, :math:`T` (K).
          )pbdoc")
      .def_readonly("beta", &conditions_type::beta, R"pbdoc(
          float: Equal to :math:`1/(k_B T)`.
          )pbdoc")
      .def(
          "set_temperature",
          [](conditions_type &self, double temperature) {
            self.set_temperature(temperature);
          },
          R"pbdoc(
          Set temperature and update `beta` to be consistent

          Parameters
          ----------
          temperature: float
              The temperature, :math:`T` (K).
          )pbdoc",
          py::arg("temperature"))
      .def(
          "set_temperature_from_value_map",
          [](conditions_type &self, monte::ValueMap const &map) {
            self.set_temperature(map);
          },
          R"pbdoc(
          Set temperature and update `beta` to be consistent

          Parameters
          ----------
          value_map: libcasm.monte.ValueMap
              Conditions as a :class:`~libcasm.monte.ValueMap` that contains
              the scalar value ``"temperature"`` in K.
          )pbdoc",
          py::arg("value_map"))
      .def("put_temperature", &conditions_type::put_temperature,
           R"pbdoc(
          Put the temperature in a conditions ValueMap

          Parameters
          ----------
          value_map: libcasm.monte.ValueMap
              Conditions as a :class:`~libcasm.monte.ValueMap` that is updated
              to contain the scalar value ``"temperature"`` in K.
          )pbdoc",
           py::arg("value_map"))
      .def_readonly("composition_converter",
                    &conditions_type::composition_converter, R"pbdoc(
          libcasm.composition.CompositionConverter: Converter between number \
          of species per unit cell and parametric composition.
          )pbdoc")
      .def_readonly("param_chem_pot", &conditions_type::param_chem_pot, R"pbdoc(
          np.ndarray: The parametric composition, :math:`\vec{x}`.
          )pbdoc")
      .def_readonly("exchange_chem_pot", &conditions_type::exchange_chem_pot,
                    R"pbdoc(
          np.ndarray: 2d array, ``M``, with values ``M[i,j] = mu_i - mu_j``.

          The exchange chemical potential for two species, `i` and `j`, which
          are indices into `composition_converter.components()`.
          change in potential per supercell due to a change of species, \
          using indices into :py:attr:`~libcasm.clexmonte.System.species_list`,
          according to:

          .. code-block:: Python

              delta_potential_energy -= \
                  conditions.exchange_chem_pot[new_species_index, curr_species_index]


          )pbdoc")
      .def(
          "set_param_chem_pot",
          [](conditions_type &self, Eigen::VectorXd const &param_chem_pot) {
            self.set_param_chem_pot(param_chem_pot);
          },
          R"pbdoc(
          Set the parametric chemical potential and update `exchange_chem_pot`\
          to be consistent

          Parameters
          ----------
          param_chem_pot: np.array
              The parametric composition, :math:`\vec{x}`.
          )pbdoc",
          py::arg("param_chem_pot"))
      .def(
          "set_param_chem_pot_from_value_map",
          [](conditions_type &self, monte::ValueMap const &map) {
            self.set_param_chem_pot(map);
          },
          R"pbdoc(
          Set the parametric chemical potential and update `exchange_chem_pot`\
          to be consistent

          Parameters
          ----------
          param_chem_pot: np.array
              Conditions as a :class:`~libcasm.monte.ValueMap` that contains
              the vector value ``"param_chem_pot"`` with the parametric
              composition, :math:`\vec{x}`.
          )pbdoc",
          py::arg("value_map"))
      .def("put_param_chem_pot", &conditions_type::put_param_chem_pot,
           R"pbdoc(
          Put the parametric chemical potential in a conditions ValueMap

          Parameters
          ----------
          value_map: libcasm.monte.ValueMap
              Conditions as a :class:`~libcasm.monte.ValueMap` that is updated
              to contain the vector value ``"param_chem_pot"``.
          )pbdoc",
           py::arg("value_map"))
      .def("set_all", &conditions_type::set_all, py::arg("value_map"),
           py::arg("is_increment"))
      .def("to_value_map", &conditions_type::to_value_map,
           py::arg("is_increment"));

  // TODO: Document methods
  py::class_<potential_type, std::shared_ptr<potential_type>>(
      m, "SemiGrandCanonicalPotential",
      R"pbdoc(
      Semi-grand canonical potential calculator

      )pbdoc")
      .def(py::init<std::shared_ptr<system_type>>(),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        system : System
            System data.
        )pbdoc",
           py::arg("system"))
      .def("set", &potential_type::set, R"pbdoc(
          Set the current state being calculated

          Parameters
          ----------
          state: libcasm.clexmonte.MonteCarloState
              The state to be calculated
          conditions: libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalConditions
              The semi-grand canonical thermodynamic conditions. They should
              be equivalent to `state.conditions`, but converted to
              :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalConditions`
              for efficient evaluation.

          )pbdoc",
           py::arg("state"), py::arg("conditions"))
      .def_property_readonly("state", &potential_type::state, R"pbdoc(
          Optional[libcasm.clexmonte.MonteCarloState]: The current state being \
          calculated
          )pbdoc")
      .def_property_readonly("conditions", &potential_type::conditions, R"pbdoc(
          libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalConditions: The \
          semi-grand canonical thermodynamic conditions
          )pbdoc")
      .def_property_readonly("formation_energy",
                             &potential_type::formation_energy, R"pbdoc(
          libcasm.clexulator.ClusterExpansion: The formation energy cluster \
          expansion calculator
          )pbdoc")
      .def("per_supercell", &potential_type::per_supercell, R"pbdoc(
          Calculate and return the potential per supercell

          Returns
          -------
          value: float
              The potential per supercell for the current state
          )pbdoc")
      .def("per_unitcell", &potential_type::per_unitcell, R"pbdoc(
          Calculate and return the potential per unit cell

          Returns
          -------
          value: float
              The potential per unit cell for the current state
          )pbdoc")
      .def("occ_delta_per_supercell", &potential_type::occ_delta_per_supercell,
           R"pbdoc(
          Calculate and return the change in potential per supercell

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

  py::class_<calculator_type, std::shared_ptr<calculator_type>>(
      m, "SemiGrandCanonicalCalculator",
      R"pbdoc(
      Implements semi-grand canonical Monte Carlo calculations
      )pbdoc")
      .def(py::init<>(&make_calculator),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        system : libcasm.clexmonte.System
            Cluster expansion model system data.
        )pbdoc",
           py::arg("system"))
      .def(
          "make_default_sampling_fixture_params",
          [](std::shared_ptr<calculator_type> &self, std::string label,
             bool write_results, bool write_trajectory, bool write_observations,
             bool write_status, std::optional<std::string> output_dir,
             std::optional<std::string> log_file,
             double log_frequency_in_s) -> sampling_fixture_params_type {
            monte::SamplingParams sampling_params;
            {
              auto &s = sampling_params;
              s.sampler_names = {"formation_energy", "potential_energy",
                                 "mol_composition", "param_composition"};
              std::string prefix;
              prefix = "order_parameter_";
              for (auto const &pair : self->system->dof_spaces) {
                s.sampler_names.push_back(prefix + pair.first);
              }
              prefix = "subspace_order_parameter_";
              for (auto const &pair : self->system->dof_subspaces) {
                s.sampler_names.push_back(prefix + pair.first);
              }
              if (write_trajectory) {
                s.do_sample_trajectory = true;
              }
            }

            monte::CompletionCheckParams<statistics_type>
                completion_check_params;
            {
              auto &c = completion_check_params;
              c.equilibration_check_f = monte::default_equilibration_check;
              c.calc_statistics_f =
                  monte::default_statistics_calculator<statistics_type>();

              converge(self->sampling_functions, completion_check_params)
                  .set_abs_precision("formation_energy", 0.001)
                  .set_abs_precision("param_composition", 0.001);
            }

            std::vector<std::string> analysis_names = {
                "heat_capacity", "mol_susc", "param_susc",
                "mol_thermochem_susc", "param_thermochem_susc"};

            return make_sampling_fixture_params(
                label, self->sampling_functions, self->json_sampling_functions,
                self->analysis_functions, sampling_params,
                completion_check_params, analysis_names, write_results,
                write_trajectory, write_observations, write_status, output_dir,
                log_file, log_frequency_in_s);
          },
          R"pbdoc(
          Construct default sampling fixture parameters

          Notes
          -----

          By default:

          - Sampling occurs linearly, by pass, with period 1, for:

            - "formation_energy",
            - "mol_composition",
            - "param_composition",
            - "potential_energy",
            - "order_parameter_<key>" (for all DoFSpace key), and
            - "subspace_order_parameter_<key>" (for all DoFSpace key).

          - Analysis functions are evaluated for:

            - "heat_capacity",
            - "mol_susc",
            - "param_susc",
            - "mol_thermochem_susc", and
            - "param_thermochem_susc".

          - Convergence of "potential_enthalpy" is set to an
            absolute precision of 0.001, and "param_composition" to 0.001.
          - Completion is checked every 100 samples, starting with the 100-th.
          - No cutoffs are set.

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
                "libcasm.clexmonte.semigrand_canonical."
                "SemiGrandCanonicalCalculator.sampling_fixture_params_from_"
                "dict"};
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
      .def(
          "run_fixture",
          [](calculator_type &self, state_type &state,
             sampling_fixture_params_type sampling_fixture_params,
             monte::OccLocation *occ_location,
             std::optional<std::shared_ptr<engine_type>> engine)
              -> std::shared_ptr<sampling_fixture_type> {
            // construct RunManager
            if (!engine.has_value()) {
              engine = std::make_shared<engine_type>();
            }
            std::vector<sampling_fixture_params_type> _params_list;
            _params_list.push_back(sampling_fixture_params);
            bool global_cutoff = true;
            auto run_manager = std::make_shared<run_manager_type>(
                engine.value(), _params_list, global_cutoff);

            // construct OccLocation if necessary
            std::unique_ptr<monte::OccLocation> tmp;
            if (!occ_location) {
              monte::Conversions const &convert =
                  get_index_conversions(*self.system, state);
              monte::OccCandidateList const &occ_candidate_list =
                  get_occ_candidate_list(*self.system, state);
              tmp = std::make_unique<monte::OccLocation>(
                  convert, occ_candidate_list, self.update_species);
              tmp->initialize(get_occupation(state));
              occ_location = tmp.get();
            }

            // run
            self.run(state, *occ_location, *run_manager);
            return run_manager->sampling_fixtures[0];
          },
          R"pbdoc(
          Perform a single run, with a single sampling fixture, evolving the \
          input state

          Parameters
          ----------
          state : libcasm.clexmonte.MonteCarloState
              The input state.
          sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
              Specifies sampling and convergence criteria and collects results
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state` and it is used and updated during the run.
              If None, a occupant location list is generated for the run.
          engine: Optional[libcasm.monte.RandomNumberEngine] = None
                Optional random number engine to use.

          Returns
          -------
          sampling_fixture: libcasm.clexmonte.SamplingFixture
              A SamplingFixture with collected results.
          )pbdoc",
          py::arg("state"), py::arg("sampling_fixture_params"),
          py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr),
          py::arg("engine") = std::nullopt)
      .def(
          "run",
          [](calculator_type &self, state_type &state,
             std::shared_ptr<run_manager_type> run_manager,
             monte::OccLocation *occ_location)
              -> std::shared_ptr<run_manager_type> {
            // construct OccLocation if necessary
            std::unique_ptr<monte::OccLocation> tmp;
            if (!occ_location) {
              monte::Conversions const &convert =
                  get_index_conversions(*self.system, state);
              monte::OccCandidateList const &occ_candidate_list =
                  get_occ_candidate_list(*self.system, state);
              tmp = std::make_unique<monte::OccLocation>(
                  convert, occ_candidate_list, self.update_species);
              tmp->initialize(get_occupation(state));
              occ_location = tmp.get();
            }

            // run
            self.run(state, *occ_location, *run_manager);
            return run_manager;
          },
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
      .def(
          "standard_sampling_functions",
          [](std::shared_ptr<calculator_type> const &self) {
            return calculator_type::standard_sampling_functions(self);
          },
          R"pbdoc(
          Construct standard semi-grand canonical sampling functions

          Returns
          -------
          sampling_functions: libcasm.monte.StateSamplingFunctionMap
              Standard state sampling functions for semi-grand canonical Monte
              Carlo calculations.
          )pbdoc")
      .def(
          "standard_json_sampling_functions",
          [](std::shared_ptr<calculator_type> const &self) {
            return calculator_type::standard_json_sampling_functions(self);
          },
          R"pbdoc(
          Construct standard semi-grand canonical JSON sampling functions

          Returns
          -------
          json_sampling_functions: libcasm.monte.jsonStateSamplingFunctionMap
              Standard JSON state sampling functions for semi-grand canonical
              Monte Carlo calculations.
          )pbdoc")
      .def(
          "standard_analysis_functions",
          [](std::shared_ptr<calculator_type> const &self) {
            return calculator_type::standard_analysis_functions(self);
          },
          R"pbdoc(
          Construct standard semi-grand canonical results analysis functions

          Returns
          -------
          analysis_functions: libcasm.clexmonte.ResultsAnalysisFunctionMap
              Standard results analysis functions for semi-grand canonical
              Monte Carlo calculations.
          )pbdoc")
      .def_readonly("system", &calculator_type::system, R"pbdoc(
          System : System data.
          )pbdoc")
      .def_readonly("update_species", &calculator_type::update_species,
                    R"pbdoc(
          bool : True if this type of calculation tracks species location \
          changes; False otherwise.
          )pbdoc")
      .def_readonly("time_sampling_allowed",
                    &calculator_type::time_sampling_allowed, R"pbdoc(
          bool : True if this calculation allows time-based sampling; \
          False otherwise.
          )pbdoc")
      .def_readonly("state", &calculator_type::state, R"pbdoc(
          Optional[MonteCarloState] : The current state.
          )pbdoc")
      .def_readonly("transformation_matrix_to_super",
                    &calculator_type::transformation_matrix_to_super,
                    R"pbdoc(
          np.ndarray[np.int64] : The current state's supercell transformation \
          matrix.
          )pbdoc")
      .def_readonly("occ_location", &calculator_type::occ_location, R"pbdoc(
          libcasm.monte.events.OccLocation : The current state's occupant \
          location list.
          )pbdoc")
      .def_readonly("conditions", &calculator_type::conditions, R"pbdoc(
          SemiGrandCanonicalConditions : The current state's conditions.
          )pbdoc")
      .def_readonly("potential", &calculator_type::potential, R"pbdoc(
          SemiGrandCanonicalPotential : The current state's potential \
          calculator.
          )pbdoc")
      .def_readonly("formation_energy", &calculator_type::formation_energy,
                    R"pbdoc(
          libcasm.clexulator.ClusterExpansion : The current state's formation \
          energy cluster expansion calculator.
          )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
