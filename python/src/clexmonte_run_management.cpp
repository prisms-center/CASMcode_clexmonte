#include <pybind11/eigen.h>
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
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/Configuration_json_io.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/monte/run_management/io/json/jsonResultsIO_impl.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// used for libcasm.clexmonte:
typedef std::mt19937_64 engine_type;
typedef monte::RandomNumberGenerator<engine_type> generator_type;
typedef clexmonte::config_type config_type;
typedef clexmonte::state_type state_type;
typedef clexmonte::statistics_type statistics_type;

typedef monte::SamplingFixture<config_type, statistics_type, engine_type>
    sampling_fixture_type;
typedef clexmonte::sampling_fixture_params_type sampling_fixture_params_type;
typedef clexmonte::run_manager_type<engine_type> run_manager_type;

typedef monte::Results<config_type, statistics_type> results_type;
typedef clexmonte::results_io_type results_io_type;
typedef monte::jsonResultsIO<results_type> json_results_io_type;

typedef monte::ResultsAnalysisFunction<config_type, statistics_type>
    analysis_function_type;
typedef monte::ResultsAnalysisFunctionMap<config_type, statistics_type>
    analysis_function_map_type;

analysis_function_type make_analysis_function(
    std::string name, std::string description, std::vector<Index> shape,
    std::function<Eigen::VectorXd(results_type const &)> function,
    std::optional<std::vector<std::string>> component_names) {
  if (function == nullptr) {
    throw std::runtime_error(
        "Error constructing ResultsAnalysisFunction: function == nullptr");
  }
  if (!component_names.has_value()) {
    return analysis_function_type(name, description, shape, function);
  } else {
    return analysis_function_type(name, description, *component_names, shape,
                                  function);
  }
}

sampling_fixture_params_type make_sampling_fixture_params(
    std::string label, monte::StateSamplingFunctionMap sampling_functions,
    monte::jsonStateSamplingFunctionMap json_sampling_functions,
    analysis_function_map_type analysis_functions,
    monte::SamplingParams sampling_params,
    monte::CompletionCheckParams<statistics_type> completion_check_params,
    std::optional<std::string> output_dir, bool write_trajectory,
    bool write_observations, std::optional<std::string> log_file,
    double log_frequency_in_s) {
  if (!output_dir.has_value()) {
    output_dir = (fs::path("output") / label).string();
  }
  if (!log_file.has_value()) {
    log_file = (fs::path(*output_dir) / "status.json").string();
  }

  std::unique_ptr<results_io_type> results_io =
      std::make_unique<json_results_io_type>(*output_dir, write_trajectory,
                                             write_observations);
  monte::MethodLog method_log;
  method_log.logfile_path = *log_file;
  method_log.log_frequency = log_frequency_in_s;
  return sampling_fixture_params_type(
      label, sampling_functions, json_sampling_functions, analysis_functions,
      sampling_params, completion_check_params, std::move(results_io),
      method_log);
}

results_type make_results(sampling_fixture_params_type const &params) {
  return results_type(
      params.sampling_params.sampler_names, params.sampling_functions,
      params.sampling_params.json_sampler_names, params.json_sampling_functions,
      params.analysis_functions);
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MAKE_OPAQUE(CASM::monte::SamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonSamplerMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::StateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASM::monte::jsonStateSamplingFunctionMap);
PYBIND11_MAKE_OPAQUE(CASMpy::analysis_function_map_type);

PYBIND11_MODULE(_clexmonte_run_management, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Cluster expansion Monte Carlo classes and methods

        libcasm.clexmonte._clexmonte_run_management
        -------------------------------------------

        The RunManager class:

        - holds sampling fixtures,
        - checks for completion, and
        - collects results.
    )pbdoc";
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.sampling");

  py::class_<results_type> pyResults(m, "Results",
                                     R"pbdoc(
        Data structure that collects Monte Carlo results from one sampling \
        fixture.
        )pbdoc");

  py::class_<analysis_function_type>(m, "ResultsAnalysisFunction",
                                     R"pbdoc(
        Calculates functions of the sampled data at the end of a run

        )pbdoc")
      .def(py::init<>(&make_analysis_function),
           R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          name : str
              Name of the sampled quantity.
          description : str
              Description of the function.
          shape : List[int]
              Shape of quantity, with column-major unrolling

              Scalar: [], Vector: [n], Matrix: [m, n], etc.

          function : function
              A function of :class:`~libcasm.clexmonte.Results` that returns
              an array of the proper size.
          component_names : Optional[List[str]] = None
              A name for each component of the resulting vector.

              Can be strings representing an indices (i.e "0", "1", "2", etc.)
              or can be a descriptive string (i.e. "Mg", "Va", "O", etc.). If
              None, indices for column-major ordering are used (i.e. "0,0",
              "1,0", ..., "m-1,n-1")

          )pbdoc",
           py::arg("name"), py::arg("description"), py::arg("shape"),
           py::arg("function"), py::arg("component_names") = std::nullopt)
      .def_readwrite("name", &analysis_function_type::name,
                     R"pbdoc(
          str : Name of the analysis function.
          )pbdoc")
      .def_readwrite("description", &analysis_function_type::description,
                     R"pbdoc(
          str : Description of the function.
          )pbdoc")
      .def_readwrite("shape", &analysis_function_type::shape,
                     R"pbdoc(
          List[int] : Shape of quantity, with column-major unrolling.

          Scalar: [], Vector: [n], Matrix: [m, n], etc.
          )pbdoc")
      .def_readwrite("component_names",
                     &analysis_function_type::component_names,
                     R"pbdoc(
          List[str] : A name for each component of the resulting vector.

          Can be strings representing an indices (i.e "0", "1", "2", etc.) or can be a descriptive string (i.e. "Mg", "Va", "O", etc.). If the sampled quantity is an unrolled matrix, indices for column-major ordering are typical (i.e. "0,0", "1,0", ..., "m-1,n-1").
          )pbdoc")
      .def_readwrite("function", &analysis_function_type::function,
                     R"pbdoc(
          function : The function to be evaluated.

          A function of :class:`~libcasm.clexmonte.Results` that returns
          an array of the proper size.
          )pbdoc")
      .def(
          "__call__",
          [](analysis_function_type const &f, results_type const &results)
              -> Eigen::VectorXd { return f(results); },
          R"pbdoc(
          Evaluates the function

          Equivalent to calling :py::attr:`~libcasm.monte.ResultsAnalysisFunction.function`.
          )pbdoc");

  py::bind_map<analysis_function_map_type>(m, "ResultsAnalysisFunctionMap",
                                           R"pbdoc(
    ResultsAnalysisFunctionMaP stores :class:`~libcasm.clexmonte.ResultsAnalysisFunction` by name.

    Notes
    -----
    ResultsAnalysisFunctionMap is a Dict[str, :class:`~libcasm.clexmonteResultsAnalysisFunction`]-like object.
    )pbdoc",
                                           py::module_local(false));

  py::class_<sampling_fixture_params_type>(m, "SamplingFixtureParams",
                                           R"pbdoc(
      Sampling fixture parameters

      Specifies what to sample, when, and how to check for completion.
      )pbdoc")
      .def(py::init<>(&make_sampling_fixture_params),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          label: str
              Label for the :class:`SamplingFixture`.
          sampling_functions : libcasm.monte.sampling.StateSamplingFunctionMap
              All possible state sampling functions
          json_sampling_functions: libcasm.monte.sampling.StateSamplingFunctionMap
              All possible JSON state sampling functions
          analysis_functions: ResultsAnalysisFunctionMap
              Results analysis functions
          sampling_params: libcasm.monte.sampling.SamplingParams
              Sampling parameters, specifies which sampling functions to call
          completion_check_params: libcasm.monte.sampling.CompletionCheckParams
              Completion check parameters
          output_dir: Optional[str] = None
              Directory in which write results. If None, uses
              ``"output" / label``.
          write_trajectory: bool = False
              If True, write the trajectory of Monte Carlo states when each
              sample taken.
          write_observations: bool = False
              If True, write all individual sample observations. Otherwise, only
              mean and estimated precision are written.
          log_file: str = Optional[str] = None
              Path to where a run status log file should be written with run
              information. If None, uses ``output_dir / "status.json"``.
          log_frequency_in_s: float = 600.0
              Minimum time between when the status log should be written, in
              seconds. The status log is only written after a sample is taken,
              so if the `sampling_params` are such that the time between
              samples is longer than `log_frequency_is_s` the status log will
              be written less frequently.
          )pbdoc",
           py::arg("label"), py::arg("sampling_functions"),
           py::arg("json_sampling_functions"), py::arg("analysis_functions"),
           py::arg("sampling_params"), py::arg("completion_check_params"),
           py::arg("output_dir") = std::nullopt,
           py::arg("write_trajectory") = false,
           py::arg("write_observations") = false,
           py::arg("log_file") = std::nullopt,
           py::arg("log_frequency_in_s") = 600.0)
      .def_static(
          "from_dict",
          [](const nlohmann::json &data, std::string label,
             monte::StateSamplingFunctionMap const &sampling_functions,
             monte::jsonStateSamplingFunctionMap const &json_sampling_functions,
             analysis_function_map_type const &analysis_functions,
             bool time_sampling_allowed) {
            jsonParser json{data};
            InputParser<sampling_fixture_params_type> parser(
                json, label, sampling_functions, json_sampling_functions,
                analysis_functions, clexmonte::standard_results_io_methods(),
                time_sampling_allowed);
            std::runtime_error error_if_invalid{
                "Error in libcasm.clexmonte.SamplingFixtureParams.from_dict"};
            report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
            return std::move(*parser.value);
          },
          R"pbdoc(
          Construct SamplingFixtureParams from a Python dict

          Parameters
          ----------
          data: dict
              The input data
          label: str
              Label for the :class:`SamplingFixture`.
          sampling_functions : libcasm.monte.sampling.StateSamplingFunctionMap
              All possible state sampling functions
          json_sampling_functions: libcasm.monte.sampling.StateSamplingFunctionMap
              All possible JSON state sampling functions
          analysis_functions: ResultsAnalysisFunctionMap
              Results analysis functions
          time_sampling_allowed: bool
              Validates input based on whether the intended Monte Carlo
              calculator allows time-based sampling or not.

          Returns
          -------
          sampling_fixture_params: SamplingFixture
              The SamplingFixtureParams
          )pbdoc",
          py::arg("data"), py::arg("label"), py::arg("sampling_functions"),
          py::arg("json_sampling_functions"), py::arg("analysis_functions"),
          py::arg("time_sampling_allowed"));

  pyResults
      .def(py::init<>(&make_results),
           R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          sampling_fixture_params: SamplingFixtureParams
              Sampling fixture parameters.
          )pbdoc",
           py::arg("sampling_fixture_params"))
      .def_readwrite("sampler_names", &results_type::sampler_names,
                     R"pbdoc(
          list[str] : The names of sampling functions that will be sampled. \
              Must be keys in `sampling_functions`.
          )pbdoc")
      .def_readwrite("sampling_functions", &results_type::sampling_functions,
                     R"pbdoc(
          libcasm.monte.StateSamplingFunctionMap : State sampling functions.
          )pbdoc")
      .def_readwrite("json_sampler_names", &results_type::json_sampler_names,
                     R"pbdoc(
          list[str] : The names of JSON sampling functions that will be \
              sampled. Must be keys in `json_sampling_functions`.
          )pbdoc")
      .def_readwrite("json_sampling_functions",
                     &results_type::json_sampling_functions,
                     R"pbdoc(
          libcasm.monte.jsonStateSamplingFunctionMap : JSON state sampling functions.
          )pbdoc")
      .def_readwrite("analysis_functions", &results_type::analysis_functions,
                     R"pbdoc(
          libcasm.clexmonte.ResultsAnalysisFunctionMap : Results analysis \
          functions. All will be evaluated.
          )pbdoc")
      .def_readwrite("elapsed_clocktime", &results_type::elapsed_clocktime,
                     R"pbdoc(
          Optional[float] : Elapsed clocktime
          )pbdoc")
      .def_readwrite("samplers", &results_type::samplers,
                     R"pbdoc(
          libcasm.monte.sampling.SamplerMap : Sampled data
          )pbdoc")
      .def_readwrite("json_samplers", &results_type::json_samplers,
                     R"pbdoc(
          libcasm.monte.sampling.jsonSamplerMap : JSON sampled data
          )pbdoc")
      .def_readwrite("analysis", &results_type::analysis,
                     R"pbdoc(
          dict[str, np.ndarray] : Results of analysis functions
          )pbdoc")
      .def_readwrite("sample_count", &results_type::sample_count,
                     R"pbdoc(
          list[int] : Count (could be number of passes or number of steps) when\
          samples occurred
          )pbdoc")
      .def_readwrite("sample_time", &results_type::sample_time,
                     R"pbdoc(
          list[float] : Time when samples occurred (if applicable, may be \
          empty).
          )pbdoc")
      .def_readwrite("sample_weight", &results_type::sample_weight,
                     R"pbdoc(
          list[float] : Weights given to samples (not normalized, may be empty).
          )pbdoc")
      .def_readwrite("sample_clocktime", &results_type::sample_clocktime,
                     R"pbdoc(
          list[float] : Elapsed clocktime when a sample occurred.
          )pbdoc")
      .def_readwrite("sample_trajectory", &results_type::sample_trajectory,
                     R"pbdoc(
          list[MonteCarloConfiguration] : Configuration when a sample occurred \
          (if requested, may be empty).
          )pbdoc")
      .def_readwrite("completion_check_results",
                     &results_type::completion_check_results,
                     R"pbdoc(
          list[libcasm.monte.CompletionCheckParams] : Completion check results
          )pbdoc")
      .def_readwrite("n_accept", &results_type::n_accept,
                     R"pbdoc(
          int : Total number of acceptances
          )pbdoc")
      .def_readwrite("n_reject", &results_type::n_reject,
                     R"pbdoc(
          int : Total number of rejections
          )pbdoc");

  py::class_<sampling_fixture_type>(m, "SamplingFixture",
                                    R"pbdoc(
      Sampling fixture

      A data structure that collects sampled data during a Monte Carlo run and
      completion check results.
      )pbdoc")
      .def(py::init<sampling_fixture_params_type const &,
                    std::shared_ptr<engine_type>>(),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          sampling_fixture_params: list[SamplingFixtureParams]
              Sampling fixture parameters, specifying what to sample, when, and
              how to check for completion.
          engine : libcasm.monte.RandomNumberEngine
              Random number generation engine
          )pbdoc",
           py::arg("sampling_fixture_params"), py::arg("engine"));

  py::class_<run_manager_type>(m, "RunManager",
                               R"pbdoc(
      RunManager is a collection one or more SamplingFixture given to a \
      Monte Carlo run method.

      )pbdoc")
      .def(py::init<std::shared_ptr<engine_type>,
                    std::vector<sampling_fixture_params_type>, bool>(),
           R"pbdoc(
          .. rubric:: Constructor

          Parameters
          ----------
          engine : libcasm.monte.RandomNumberEngine
              Random number generation engine
          sampling_fixture_params: list[SamplingFixtureParams]
              Sampling fixture parameters, specifying what to sample, when, and
              how to check for completion.
          global_cutoff: bool = True
              If true, the run is complete if any sampling fixture is complete.
              Otherwise, all sampling fixtures must be completed for the run
              to be completed.
          )pbdoc",
           py::arg("engine"), py::arg("sampling_fixture_params"),
           py::arg("global_cutoff") = true);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
