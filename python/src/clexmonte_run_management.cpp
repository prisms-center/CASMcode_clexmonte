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
#include "casm/clexmonte/run/functions.hh"
#include "casm/clexmonte/run/io/json/RunParams_json_io.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/monte/run_management/RunManager.hh"
#include "casm/monte/run_management/io/json/SamplingFixtureParams_json_io.hh"
#include "casm/monte/run_management/io/json/jsonResultsIO_impl.hh"
#include "casm/monte/sampling/Sampler.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// used for libcasm.clexmonte:
typedef clexmonte::default_engine_type engine_type;
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

results_type make_results(sampling_fixture_params_type const &params) {
  return results_type(
      params.sampling_params.sampler_names, params.sampling_functions,
      params.sampling_params.json_sampler_names, params.json_sampling_functions,
      params.analysis_functions);
}

monte::SAMPLE_METHOD sample_method_from_string(std::string method) {
  if (method == "linear") {
    return monte::SAMPLE_METHOD::LINEAR;
  } else if (method == "log") {
    return monte::SAMPLE_METHOD::LOG;
  } else if (method == "custom") {
    return monte::SAMPLE_METHOD::CUSTOM;
  } else {
    throw std::runtime_error(
        "Error in sample_method_from_string: invalid method '" + method + "'");
  }
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
    Cluster expansion Monte Carlo run management
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

          Equivalent to calling :py::attr:`~libcasm.clexmonte.ResultsAnalysisFunction.function`.
          )pbdoc");

  py::bind_map<analysis_function_map_type>(m, "ResultsAnalysisFunctionMap",
                                           R"pbdoc(
    ResultsAnalysisFunctionMaP stores :class:`~libcasm.clexmonte.ResultsAnalysisFunction` by name.

    Notes
    -----
    ResultsAnalysisFunctionMap is a Dict[str, :class:`~libcasm.clexmonte.ResultsAnalysisFunction`]-like object.
    )pbdoc",
                                           py::module_local(false));

  py::class_<sampling_fixture_params_type>(m, "SamplingFixtureParams",
                                           R"pbdoc(
      Sampling fixture parameters

      Specifies what to sample, when, and how to check for completion.
      )pbdoc")
      .def(py::init<>(&clexmonte::make_sampling_fixture_params),
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
          analysis_names: list[str] = []
              List of which analysis functions should be evaluated.
          write_results: bool = True
              If True, write results to file upon completion.
          write_trajectory: bool = False
              If True, write the trajectory of Monte Carlo states when each
              sample taken.
          write_observations: bool = False
              If True, write all individual sample observations. Otherwise, only
              mean and estimated precision are written.
          write_status: bool = True
              If True, write log files with convergence status.
          output_dir: Optional[str] = None
              Directory in which write results. If None, uses
              ``"output" / label``.
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
           py::arg("analysis_names") = std::vector<std::string>(),
           py::arg("write_results") = true, py::arg("write_trajectory") = false,
           py::arg("write_observations") = false,
           py::arg("write_status") = true, py::arg("output_dir") = std::nullopt,
           py::arg("log_file") = std::nullopt,
           py::arg("log_frequency_in_s") = 600.0)
      .def_readwrite("label", &sampling_fixture_params_type::label, R"pbdoc(
          str: Label, to name output and distinguish multiple sampling fixtures
          )pbdoc")
      .def_readwrite("sampling_functions",
                     &sampling_fixture_params_type::sampling_functions,
                     R"pbdoc(
          libcasm.monte.sampling.StateSamplingFunctionMap: State sampling \
          functions
          )pbdoc")
      .def_readwrite("json_sampling_functions",
                     &sampling_fixture_params_type::json_sampling_functions,
                     R"pbdoc(
          libcasm.monte.sampling.jsonStateSamplingFunctionMap: JSON state \
          sampling functions
          )pbdoc")
      .def_readwrite("analysis_functions",
                     &sampling_fixture_params_type::analysis_functions,
                     R"pbdoc(
          ResultsAnalysisFunctionMap: Results analysis functions
          )pbdoc")
      .def_readwrite("analysis_names",
                     &sampling_fixture_params_type::analysis_functions,
                     R"pbdoc(
          list[str]: List of which analysis functions should be evaluated.

          Note that this is a property that either (i) gets a copy of the
          list, or (ii) sets the entire list. Doing `x.analysis_names.append(y)`
          or `x.analysis_names += [y, z]` will not modify the SamplingFixtureParams
          object. Instead, use `x.append_to_analysis_names(y)`,
          `x.remove_from_analysis_names(y)`, or `x.extend_analysis_names([y, z])`.
          )pbdoc")
      .def(
          "append_to_analysis_names",
          [](sampling_fixture_params_type &self, std::string name) {
            self.analysis_names.push_back(name);
          },
          R"pbdoc(
          Append a name to `analysis_names`.

          Parameters
          ----------
          name: str
              Name to append to `analysis_names`.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "remove_from_analysis_names",
          [](sampling_fixture_params_type &self, std::string name) {
            if (auto it = std::find(self.analysis_names.begin(),
                                    self.analysis_names.end(), name);
                it != self.analysis_names.end()) {
              self.analysis_names.erase(it);
            }
            return self;
          },
          R"pbdoc(
          Remove a name from `analysis_names`.

          Parameters
          ----------
          name: str
              Name to remove from `analysis_names`.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "extend_analysis_names",
          [](sampling_fixture_params_type &self,
             std::vector<std::string> names) {
            self.analysis_names.insert(self.analysis_names.end(), names.begin(),
                                       names.end());
          },
          R"pbdoc(
          Append multiple names to `analysis_names`.

          Parameters
          ----------
          names: list[str]
              Names to append to `analysis_names`.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("names"))
      .def_readwrite("sampling_params",
                     &sampling_fixture_params_type::sampling_params,
                     R"pbdoc(
          libcasm.monte.sampling.SamplingParams: Sampling parameters
          )pbdoc")
      .def_readwrite("completion_check_params",
                     &sampling_fixture_params_type::completion_check_params,
                     R"pbdoc(
          libcasm.monte.sampling.CompletionCheckParams: Completion check parameters
          )pbdoc")
      .def_readwrite("completion_check_params",
                     &sampling_fixture_params_type::completion_check_params,
                     R"pbdoc(
          libcasm.monte.sampling.CompletionCheckParams: Completion check parameters
          )pbdoc")
      .def_readwrite("completion_check_params",
                     &sampling_fixture_params_type::completion_check_params,
                     R"pbdoc(
          libcasm.monte.sampling.CompletionCheckParams: Completion check parameters
          )pbdoc")
      .def(
          "set_io_params",
          [](sampling_fixture_params_type &self, bool write_results,
             bool write_trajectory, bool write_observations, bool write_status,
             std::optional<std::string> output_dir,
             std::optional<std::string> log_file, double log_frequency_in_s) {
            if (!output_dir.has_value()) {
              output_dir = (fs::path("output") / self.label).string();
            }
            if (!log_file.has_value()) {
              log_file = (fs::path(*output_dir) / "status.json").string();
            }

            std::unique_ptr<results_io_type> results_io;
            if (write_results) {
              results_io = std::make_unique<monte::jsonResultsIO<results_type>>(
                  *output_dir, write_trajectory, write_observations);
              if (write_trajectory) {
                self.sampling_params.do_sample_trajectory = true;
              }
            }
            self.results_io = std::move(results_io);

            monte::MethodLog method_log;
            if (write_status) {
              method_log.logfile_path = *log_file;
              method_log.log_frequency = log_frequency_in_s;
            }
            self.method_log = method_log;

            return self;
          },
          R"pbdoc(
          Set the results and status log IO parameters

          Parameters
          ----------
          write_results: bool = True
              If True, write results to file upon completion.
          write_trajectory: bool = False
              If True, write the trajectory of Monte Carlo states when each
              sample taken.
          write_observations: bool = False
              If True, write all individual sample observations. Otherwise, only
              mean and estimated precision are written.
          write_status: bool = True
              If True, write log files with convergence status.
          output_dir: Optional[str] = None
              Directory in which write results. If None, uses
              ``"output" / label``.
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
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("write_results") = true, py::arg("write_trajectory") = false,
          py::arg("write_observations") = false, py::arg("write_status") = true,
          py::arg("output_dir") = std::nullopt,
          py::arg("log_file") = std::nullopt,
          py::arg("log_frequency_in_s") = 600.0)
      .def(
          "converge",
          [](sampling_fixture_params_type &self, std::string quantity,
             std::optional<double> abs, std::optional<double> rel,
             std::optional<std::vector<std::string>> component_name,
             std::optional<std::vector<Index>> component_index) {
            auto it = self.sampling_functions.find(quantity);
            if (it == self.sampling_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.converge: '" << quantity
                  << "' is not in sampling_functions";
              throw std::runtime_error(msg.str());
            }
            auto const &f = it->second;

            if (!rel.has_value() && !abs.has_value()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.converge: No abs or rel "
                     "precision specified";
              throw std::runtime_error(msg.str());
            }

            if (!component_index.has_value()) {
              component_index = std::vector<Index>();
            }

            if (component_name.has_value()) {
              auto begin = f.component_names.begin();
              auto end = f.component_names.end();
              for (std::string n : component_name.value()) {
                auto it = std::find(begin, end, n);
                if (it == end) {
                  std::stringstream msg;
                  msg << "Error in SamplingFixtureParams.converge: '" << n
                      << "' is not a component of '" << quantity << "'";
                  throw std::runtime_error(msg.str());
                }
                component_index->push_back(std::distance(begin, it));
              }
            }

            // remove any duplicates and sort
            std::set<Index> _tmp(component_index->begin(),
                                 component_index->end());
            std::vector<Index> indices(_tmp.begin(), _tmp.end());

            // if nothing provided, default is to set all components
            if (indices.empty()) {
              for (Index i = 0; i < f.component_names.size(); ++i) {
                indices.push_back(i);
              }
            }

            for (Index i : indices) {
              monte::SamplerComponent key(quantity, i, f.component_names[i]);
              if (rel.has_value() && abs.has_value()) {
                self.completion_check_params.requested_precision[key] =
                    monte::RequestedPrecision::abs_and_rel(abs.value(),
                                                           rel.value());
              } else if (abs.has_value()) {
                self.completion_check_params.requested_precision[key] =
                    monte::RequestedPrecision::abs(abs.value());
              } else if (rel.has_value()) {
                self.completion_check_params.requested_precision[key] =
                    monte::RequestedPrecision::rel(rel.value());
              } else {
                std::stringstream msg;
                msg << "Error in SamplingFixtureParams.converge: Invalid "
                       "request for "
                       "unknown reason";
                throw std::runtime_error(msg.str());
              }
            }

            return self;
          },
          R"pbdoc(
          Set requested precision level for equilibration and convergence

          Allows setting absolute or relative precision to the specified level for
          the specified quantities. By default, all components are converged to
          the same level. If `component_name` or `component_index` are specified,
          then only the specified components are requested to converge to that level.

          Parameters
          ----------
          quantity: str
              The name of the quantity to be converged. Must match
              a state sampling function name.
          abs: Optional[float]=None
              The requested absolute convergence level
          rel: Optional[float]=None,
              The requested relative convergence level
          component_name: Optional[list[str]]=None
              The name of components to converge. Must be in the
              `component_names` of the state sampling function for
              the named quantity.
          component_index: Optional[list[int]]=None
              The indices of components to converge.


          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("quantity"), py::arg("abs") = std::nullopt,
          py::arg("rel") = std::nullopt,
          py::arg("component_name") = std::nullopt,
          py::arg("component_index") = std::nullopt)
      .def(
          "do_not_converge",
          [](sampling_fixture_params_type &self, std::string quantity,
             std::optional<std::vector<std::string>> component_name,
             std::optional<std::vector<Index>> component_index) {
            auto it = self.sampling_functions.find(quantity);
            if (it == self.sampling_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.do_not_converge: '"
                  << quantity << "' is not in sampling_functions";
              throw std::runtime_error(msg.str());
            }
            auto const &f = it->second;

            if (!component_index.has_value()) {
              component_index = std::vector<Index>();
            }

            if (component_name.has_value()) {
              auto begin = f.component_names.begin();
              auto end = f.component_names.end();
              for (std::string n : component_name.value()) {
                auto it = std::find(begin, end, n);
                if (it == end) {
                  std::stringstream msg;
                  msg << "Error in SamplingFixtureParams.do_not_converge: '"
                      << n << "' is not a component of '" << quantity << "'";
                  throw std::runtime_error(msg.str());
                }
                component_index->push_back(std::distance(begin, it));
              }
            }

            // remove any duplicates and sort
            std::set<Index> _tmp(component_index->begin(),
                                 component_index->end());
            std::vector<Index> indices(_tmp.begin(), _tmp.end());

            // if nothing provided, default is to unset all components
            if (indices.empty()) {
              for (Index i = 0; i < f.component_names.size(); ++i) {
                indices.push_back(i);
              }
            }

            for (Index i : indices) {
              monte::SamplerComponent key(quantity, i, f.component_names[i]);
              self.completion_check_params.requested_precision.erase(key);
            }

            return self;
          },
          R"pbdoc(
          Remove requested precision level for equilibration and convergence

          Allows removing the requested precision level for the specified quantities.

          Parameters
          ----------
          quantity: str
              The name of the quantity to be converged. Must match
              a state sampling function name.
          component_name: Optional[list[str]]=None
              The name of components to remove convergence criteria for.
              Must be in the `component_names` of the state sampling function for
              the named quantity.
          component_index: Optional[list[int]]=None
              The indices of components to remove convergence criteria for.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("quantity"), py::arg("component_name") = std::nullopt,
          py::arg("component_index") = std::nullopt)
      .def(
          "set_completion_check_spacing",
          [](sampling_fixture_params_type &self, bool log_spacing,
             std::optional<int> check_period, std::optional<int> check_begin,
             std::optional<float> check_base, std::optional<float> check_shift,
             std::optional<int> check_period_max) {
            // set defaults
            auto &check_params = self.completion_check_params;
            if (!log_spacing) {
              if (!check_period.has_value()) {
                check_period = 100;
              }
              check_params.check_begin = check_period.value();
              check_params.check_period = check_period.value();
            } else {
              check_params.check_begin = 0;
              check_params.check_base = 10.0;
              check_params.check_shift = 2.0;
              check_params.check_period_max = 10000;
            }

            // set user chosen values
            check_params.log_spacing = log_spacing;
            if (check_begin.has_value()) {
              check_params.check_begin = check_begin.value();
            }
            if (check_period.has_value()) {
              check_params.check_period = check_period.value();
            }
            if (check_base.has_value()) {
              check_params.check_base = check_base.value();
            }
            if (check_shift.has_value()) {
              check_params.check_shift = check_shift.value();
            }
            if (check_period_max.has_value()) {
              check_params.check_period_max = check_period_max.value();
            }
            return self;
          },
          R"pbdoc(
          Set completion check spacing parameters

          Completion checking is done at specific intervals to determine if the
          all requested quantities have converged to the desired precision. If
          completion checking occurs too frequently, it may slow down the
          simulation. If it occurs too infrequently, the simulation may run
          longer than necessary.


          Parameters
          ----------
          log_spacing: bool = False
              If True, use logarithmic spacing for completion checking; else use linear
              spacing. For linear spacing, the n-th check (n=0,1,2,...) will be taken when:

              .. code-block:: Python

                  sample = check_begin + check_period * n

              For logarithmic spacing, the n-th check will be taken when:

              .. code-block:: Python

                  sample = check_begin + round( check_base ** (n + check_shift) )

              However, if check(n) - check(n-1) > `check_period_max`, then subsequent
              checks are made every `check_period_max` samples.

              For linear spacing, the default is to check for completion after `100`,
              `200`, `300`, etc. samples are taken.

              For log spacing, the default is to check for completion after `100`,
              `1000`, `10000`, `20000`, `30000`, etc. samples are taken (note the
              effect of the default ``check_period_max=10000``).

              The default value is False, for linear spacing.
          check_period:  Optional[int] = None
              The linear completion checking period. Default is 100.
          check_begin:  Optional[int] = None
              The earliest sample to begin completion checking. Default is
              `check_period` for linear spacing and 0 for log spacing.
          check_base: Optional[float] = None
              The logarithmic completion checking base. Default is 10.
          check_shift: Optional[float] = None
              The shift for the logarithmic spacing exponent. Default is 2.
          check_period_max: Optional[int] = None
              The maximum check spacing for logarithmic check spacing. Default is 10000.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("log_spacing") = false, py::arg("check_begin") = std::nullopt,
          py::arg("check_period") = std::nullopt,
          py::arg("check_base") = std::nullopt,
          py::arg("check_shift") = std::nullopt,
          py::arg("check_period_max") = std::nullopt)
      .def(
          "clear_cutoffs",
          [](sampling_fixture_params_type &self) {
            self.completion_check_params.cutoff_params = {};
            return self;
          },
          R"pbdoc(
          Clear all completion check cutoffs

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc")
      .def(
          "set_min_count",
          [](sampling_fixture_params_type &self,
             std::optional<monte::CountType> value) {
            self.completion_check_params.cutoff_params.min_count = value;
            return self;
          },
          R"pbdoc(
          Set the `min_count` cutoff parameter

          Parameters
          ----------
          min_count: Optional[int]
              Set the minimum number of steps or passes, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("min_count"))
      .def(
          "set_min_time",
          [](sampling_fixture_params_type &self,
             std::optional<monte::TimeType> value) {
            self.completion_check_params.cutoff_params.min_time = value;
            return self;
          },
          R"pbdoc(
          Set the `min_time` cutoff parameter

          Parameters
          ----------
          min_time: Optional[float]
              Set the minimum simulated time, if applicable, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("min_time"))
      .def(
          "set_min_sample",
          [](sampling_fixture_params_type &self,
             std::optional<monte::CountType> value) {
            self.completion_check_params.cutoff_params.min_sample = value;
            return self;
          },
          R"pbdoc(
          Set the `min_sample` cutoff parameter

          Parameters
          ----------
          min_sample: Optional[int]
              Set the minimum number of samples, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("min_sample"))
      .def(
          "set_min_clocktime",
          [](sampling_fixture_params_type &self,
             std::optional<monte::TimeType> value) {
            self.completion_check_params.cutoff_params.min_clocktime = value;
            return self;
          },
          R"pbdoc(
          Set the `min_clocktime` cutoff parameter

          Parameters
          ----------
          min_clocktime: Optional[float]
              Set the minimum elapsed clocktime, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("min_clocktime"))
      .def(
          "set_max_count",
          [](sampling_fixture_params_type &self,
             std::optional<monte::CountType> value) {
            self.completion_check_params.cutoff_params.max_count = value;
            return self;
          },
          R"pbdoc(
          Set the `max_count` cutoff parameter

          Parameters
          ----------
          max_count: Optional[int]
              Set the maximum number of steps or passes, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("max_count"))
      .def(
          "set_max_time",
          [](sampling_fixture_params_type &self,
             std::optional<monte::TimeType> value) {
            self.completion_check_params.cutoff_params.max_time = value;
            return self;
          },
          R"pbdoc(
          Set the `max_time` cutoff parameter

          Parameters
          ----------
          max_time: Optional[float]
              Set the maximum simulated time, if applicable, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("max_time"))
      .def(
          "set_max_sample",
          [](sampling_fixture_params_type &self,
             std::optional<monte::CountType> value) {
            self.completion_check_params.cutoff_params.max_sample = value;
            return self;
          },
          R"pbdoc(
          Set the `max_sample` cutoff parameter

          Parameters
          ----------
          max_sample: Optional[int]
              Set the maximum number of samples, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("max_sample"))
      .def(
          "set_max_clocktime",
          [](sampling_fixture_params_type &self,
             std::optional<monte::TimeType> value) {
            self.completion_check_params.cutoff_params.max_clocktime = value;
            return self;
          },
          R"pbdoc(
          Set the `max_clocktime` cutoff parameter

          Parameters
          ----------
          max_clocktime: Optional[float]
              Set the maximum elapsed clocktime, or if None, remove the
              cutoff.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("max_clocktime"))
      .def_readwrite("method_log", &sampling_fixture_params_type::method_log,
                     R"pbdoc(
          libcasm.monte.MethodLog: Handles status file output
          )pbdoc")
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
          py::arg("time_sampling_allowed"))
      .def(
          "sample_by_step",
          [](sampling_fixture_params_type &self, std::optional<double> begin,
             double period, std::string method, double base, double shift,
             std::optional<std::function<double(monte::CountType)>>
                 custom_sample_at) {
            self.sampling_params.sample_mode = monte::SAMPLE_MODE::BY_STEP;
            self.sampling_params.sample_method =
                sample_method_from_string(method);
            if (method == "log") {
              self.sampling_params.begin = begin.value_or(0.0);
            } else {
              self.sampling_params.begin = begin.value_or(period);
            }
            self.sampling_params.period = period;
            self.sampling_params.sample_method =
                sample_method_from_string(method);
            self.sampling_params.base = base;
            self.sampling_params.shift = shift;
            if (custom_sample_at.has_value()) {
              self.sampling_params.custom_sample_at = custom_sample_at.value();
            }
            return self;
          },
          R"pbdoc(
          Set the sampling mode to "by step" and set sample spacing

          See :class:`~libcasm.monte.sampling.SamplingParams` for more details
          on how the sample spacing is calculated.

          Parameters
          ----------
          begin: Optional[float] = None
              The `begin` sample spacing parameter. If None, uses `period`
              if method != "log", and uses 0.0 if method == "log".
          period: float = 1.0
              The `period` sampling spacing parameter
          method: str = "linear"
              One of "linear", "log", or "custom"
          base: float = math.pow(10.0, 1.0/10.0)
              Base for log sampling (ignored if method != "log")
          shift: float = 1.0
              Shift for log sampling (ignored if method != "log")
          custom_sample_at: Optional[Callable[[int], float]] = None
              Custom function for sampling, which must have the signature
              ``def custom_sample_at(n: int) -> float``, returning the
              step when the n-th sample sample should be taken
              (ignored if method != "custom")

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("begin") = std::nullopt, py::arg("period") = 1.0,
          py::arg("method") = "linear",
          py::arg("base") = std::pow(10.0, 1.0 / 10.0), py::arg("shift") = 1.0,
          py::arg("custom_sample_at") = std::nullopt)
      .def(
          "sample_by_pass",
          [](sampling_fixture_params_type &self, std::optional<double> begin,
             double period, std::string method, double base, double shift,
             std::optional<std::function<double(monte::CountType)>>
                 custom_sample_at) {
            self.sampling_params.sample_mode = monte::SAMPLE_MODE::BY_PASS;
            self.sampling_params.sample_method =
                sample_method_from_string(method);
            if (method == "log") {
              self.sampling_params.begin = begin.value_or(0.0);
            } else {
              self.sampling_params.begin = begin.value_or(period);
            }
            self.sampling_params.period = period;
            self.sampling_params.base = base;
            self.sampling_params.shift = shift;
            if (custom_sample_at.has_value()) {
              self.sampling_params.custom_sample_at = custom_sample_at.value();
            }
            return self;
          },
          R"pbdoc(
          Set the sampling mode to "by pass" and set sample spacing

          See :class:`~libcasm.monte.sampling.SamplingParams` for more details
          on how the sample spacing is calculated.

          Parameters
          ----------
          begin: Optional[float] = None
              The `begin` sample spacing parameter. If None, uses `period`
              if method != "log", and uses 0.0 if method == "log".
          period: float = 1.0
              The `period` sampling spacing parameter
          method: str = "linear"
              One of "linear", "log", or "custom"
          base: float = math.pow(10.0, 1.0/10.0)
              Base for log sampling (ignored if method != "log")
          shift: float = 1.0
              Shift for log sampling (ignored if method != "log")
          custom_sample_at: Optional[Callable[[int], float]] = None
              Custom function for sampling, which must have the signature
              ``def custom_sample_at(n: int) -> float``, returning the
              pass when the n-th sample sample should be taken
              (ignored if method != "custom")

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("begin") = std::nullopt, py::arg("period") = 1.0,
          py::arg("method") = "linear",
          py::arg("base") = std::pow(10.0, 1.0 / 10.0), py::arg("shift") = 1.0,
          py::arg("custom_sample_at") = std::nullopt)
      .def(
          "sample_by_time",
          [](sampling_fixture_params_type &self, std::optional<double> begin,
             double period, std::string method, double base, double shift,
             std::optional<std::function<double(monte::CountType)>>
                 custom_sample_at) {
            self.sampling_params.sample_mode = monte::SAMPLE_MODE::BY_TIME;
            self.sampling_params.sample_method =
                sample_method_from_string(method);
            if (method == "log") {
              self.sampling_params.begin = begin.value_or(0.0);
            } else {
              self.sampling_params.begin = begin.value_or(period);
            }
            self.sampling_params.period = period;
            self.sampling_params.base = base;
            self.sampling_params.shift = shift;
            if (custom_sample_at.has_value()) {
              self.sampling_params.custom_sample_at = custom_sample_at.value();
            }
            return self;
          },
          R"pbdoc(
          Set the sampling mode to "by time" and set sample spacing

          See :class:`~libcasm.monte.sampling.SamplingParams` for more details
          on how the sample spacing is calculated.

          Parameters
          ----------
          begin: Optional[float] = None
              The `begin` sample spacing parameter. If None, uses `period`
              if method != "log", and uses 0.0 if method == "log".
          period: float = 1.0
              The `period` sampling spacing parameter
          method: str = "linear"
              One of "linear", "log", or "custom"
          base: float = math.pow(10.0, 1.0/10.0)
              Base for log sampling (ignored if method != "log")
          shift: float = 1.0
              Shift for log sampling (ignored if method != "log")
          custom_sample_at: Optional[Callable[[int], float]] = None
              Custom function for sampling, which must have the signature
              ``def custom_sample_at(n: int) -> float``, returning the
              pass when the n-th sample sample should be taken
              (ignored if method != "custom")

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("begin") = std::nullopt, py::arg("period") = 1.0,
          py::arg("method") = "linear",
          py::arg("base") = std::pow(10.0, 1.0 / 10.0), py::arg("shift") = 1.0,
          py::arg("custom_sample_at") = std::nullopt)
      .def(
          "sample",
          [](sampling_fixture_params_type &self, std::string quantity) {
            // If quantity is not in sampling_functions or
            // json_sampling_functions, throw an error
            if (self.sampling_functions.find(quantity) ==
                    self.sampling_functions.end() &&
                self.json_sampling_functions.find(quantity) ==
                    self.json_sampling_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.sample: '" << quantity
                  << "' is not in sampling_functions or "
                     "json_sampling_functions";
              throw std::runtime_error(msg.str());
            }
            if (self.sampling_functions.find(quantity) !=
                self.sampling_functions.end()) {
              // If not in sampling_params.sampler_names, add it
              if (std::find(self.sampling_params.sampler_names.begin(),
                            self.sampling_params.sampler_names.end(),
                            quantity) ==
                  self.sampling_params.sampler_names.end()) {
                self.sampling_params.sampler_names.push_back(quantity);
              }
            } else {
              // If not in sampling_params.json_sampler_names, add it
              if (std::find(self.sampling_params.json_sampler_names.begin(),
                            self.sampling_params.json_sampler_names.end(),
                            quantity) ==
                  self.sampling_params.json_sampler_names.end()) {
                self.sampling_params.json_sampler_names.push_back(quantity);
              }
            }
            return self;
          },
          R"pbdoc(
          Add the name of a sampling function or JSON sampling function to the
          list of quantities to be sampled

          Parameters
          ----------
          name : str
              The name of a sampling function (or JSON sampling function)
              to be added to `sampling_params.sampler_names` (or
              `sampling_params.json_sampler_names`).

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "do_not_sample",
          [](sampling_fixture_params_type &self, std::string quantity) {
            // If quantity is not in sampling_functions or
            // json_sampling_functions, throw an error
            if (self.sampling_functions.find(quantity) ==
                    self.sampling_functions.end() &&
                self.json_sampling_functions.find(quantity) ==
                    self.json_sampling_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.do_not_sample: '"
                  << quantity
                  << "' is not in sampling_functions or "
                     "json_sampling_functions";
              throw std::runtime_error(msg.str());
            }
            if (self.sampling_functions.find(quantity) !=
                self.sampling_functions.end()) {
              // If in sampling_params.sampler_names, remove it
              auto it =
                  std::find(self.sampling_params.sampler_names.begin(),
                            self.sampling_params.sampler_names.end(), quantity);
              if (it != self.sampling_params.sampler_names.end()) {
                self.sampling_params.sampler_names.erase(it);
              }
            } else {
              // If in sampling_params.json_sampler_names, remove it
              auto it = std::find(
                  self.sampling_params.json_sampler_names.begin(),
                  self.sampling_params.json_sampler_names.end(), quantity);
              if (it != self.sampling_params.json_sampler_names.end()) {
                self.sampling_params.json_sampler_names.erase(it);
              }
            }
            return self;
          },
          R"pbdoc(
          Remove the name of a sampling function or JSON sampling function from the
          list of quantities to be sampled

          Parameters
          ----------
          name : str
              The name of a sampling function (or JSON sampling function)
              to be removed from `sampling_params.sampler_names` (or
              `sampling_params.json_sampler_names`).

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "perform_analysis",
          [](sampling_fixture_params_type &self, std::string name) {
            // If name is not in analysis_functions, throw an error
            if (self.analysis_functions.find(name) ==
                self.analysis_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.perform_analysis: '"
                  << name << "' is not in analysis_functions";
              throw std::runtime_error(msg.str());
            }
            // If not in analysis_names, add it
            if (std::find(self.analysis_names.begin(),
                          self.analysis_names.end(),
                          name) == self.analysis_names.end()) {
              self.analysis_names.push_back(name);
            }
            return self;
          },
          R"pbdoc(
          Add an analysis function to the list of analyses to be performed

          Parameters
          ----------
          name : str
              The name of an analysis function included in `analysis_functions`.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "do_not_perform_analysis",
          [](sampling_fixture_params_type &self, std::string name) {
            // If name is not in analysis_functions, throw an error
            if (self.analysis_functions.find(name) ==
                self.analysis_functions.end()) {
              std::stringstream msg;
              msg << "Error in SamplingFixtureParams.do_not_perform_analysis: '"
                  << name << "' is not in analysis_functions";
              throw std::runtime_error(msg.str());
            }
            // If in analysis_names, remove it
            auto it = std::find(self.analysis_names.begin(),
                                self.analysis_names.end(), name);
            if (it != self.analysis_names.end()) {
              self.analysis_names.erase(it);
            }
            return self;
          },
          R"pbdoc(
          Remove an analysis function from the list of analyses to be performed

          Parameters
          ----------
          name : str
              The name of an analysis function to be removed from `analysis_names`.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("name"))
      .def(
          "add_sampling_function",
          [](sampling_fixture_params_type &self,
             monte::StateSamplingFunction const &f) {
            self.sampling_functions.emplace(f.name, f);
            return self;
          },
          R"pbdoc(
          Add a state sampling function

          Parameters
          ----------
          f : :class:`~libcasm.monte.sampling.StateSamplingFunction`
              The state sampling function. The function is not automatically
              added to the list of quantities to be sampled.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("f"))
      .def(
          "add_json_sampling_function",
          [](sampling_fixture_params_type &self,
             monte::jsonStateSamplingFunction const &f) {
            self.json_sampling_functions.emplace(f.name, f);
            return self;
          },
          R"pbdoc(
          Add a JSON state sampling function

          Parameters
          ----------
          f : :class:`~libcasm.monte.sampling.jsonStateSamplingFunction`
              The JSON state sampling function. The function is not automatically
              added to the list of quantities to be sampled.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("f"))
      .def(
          "add_analysis_function",
          [](sampling_fixture_params_type &self,
             analysis_function_type const &f) {
            self.analysis_functions.emplace(f.name, f);
            return self;
          },
          R"pbdoc(
          Add a results analysis function

          Parameters
          ----------
          f : :class:`~libcasm.clexmonte.ResultsAnalysisFunction`
              The results analysis function. The function is not automatically
              added to the list of analysis to be performed.

          Returns
          -------
          self: SamplingFixtureParams
              To allow chaining multiple calls, `self` is returned
          )pbdoc",
          py::arg("f"))
      .def(
          "to_dict",
          [](sampling_fixture_params_type const &self) {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the SamplingFixtureParams as a Python dict.

          Returns
          -------
          data : json
              The SamplingFixtureParams as a Python dict.
          )pbdoc");

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
          list[ libcasm.configuration.Configuration] : Configuration when a \
          sample occurred (if requested, may be empty).
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

  py::class_<sampling_fixture_type, std::shared_ptr<sampling_fixture_type>>(
      m, "SamplingFixture",
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
           py::arg("sampling_fixture_params"), py::arg("engine"))
      .def_property_readonly("label", &sampling_fixture_type::label, R"pbdoc(
          str : Label for the SamplingFixture.
          )pbdoc")
      .def_property_readonly("params", &sampling_fixture_type::params, R"pbdoc(
          SamplingFixtureParams : Access sampling fixture parameters.
          )pbdoc")
      .def_property_readonly("results", &sampling_fixture_type::results,
                             R"pbdoc(
          Results : Access sampling fixture results.
          )pbdoc")
      .def_property_readonly(
          "sample_mode",
          [](sampling_fixture_type const &self) {
            return self.counter().sample_mode;
          },
          R"pbdoc(
          libcasm.monte.sampling.SAMPLE_MODE : Sampling mode.
          )pbdoc")
      .def_property_readonly(
          "steps_per_pass",
          [](sampling_fixture_type const &self) {
            return self.counter().steps_per_pass;
          },
          R"pbdoc(
          int : Number of steps per pass (usually the number of mutating sites).
          )pbdoc")
      .def_property_readonly(
          "n_step",
          [](sampling_fixture_type const &self) { return self.counter().step; },
          R"pbdoc(
          int : Current step count (in range [0, steps_per_pass)).
          )pbdoc")
      .def_property_readonly(
          "n_pass",
          [](sampling_fixture_type const &self) { return self.counter().pass; },
          R"pbdoc(
          int : Current pass count.
          )pbdoc")
      .def_property_readonly(
          "count",
          [](sampling_fixture_type const &self) {
            return self.counter().count;
          },
          R"pbdoc(
          int : Current step or pass count, depending on `sample_mode`.
          )pbdoc");

  py::class_<run_manager_type, std::shared_ptr<run_manager_type>>(m,
                                                                  "RunManager",
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
          engine: libcasm.monte.RandomNumberEngine
              Random number generation engine to be used by
              :class:`~libcasm.clexmonte.MonteCalculator`.
          sampling_fixture_params: list[SamplingFixtureParams]
              Sampling fixture parameters, specifying what to sample, when, and
              how to check for completion.
          global_cutoff: bool = True
              If true, the run is complete if any sampling fixture is complete.
              Otherwise, all sampling fixtures must be completed for the run
              to be completed.
          )pbdoc",
           py::arg("engine"), py::arg("sampling_fixture_params"),
           py::arg("global_cutoff") = true)
      .def_readwrite("run_index", &run_manager_type::run_index,
                     R"pbdoc(
          int: Current run index
          )pbdoc")
      .def_readonly("sampling_fixtures", &run_manager_type::sampling_fixtures,
                    R"pbdoc(
          list[SamplingFixture]: List of SamplingFixture being managed.
          )pbdoc")
      .def_property_readonly(
          "sampling_fixture_labels",
          [](run_manager_type const &self) {
            std::vector<std::string> labels;
            for (auto &fixture_ptr : self.sampling_fixtures) {
              labels.push_back(fixture_ptr->label());
            }
            return labels;
          },
          R"pbdoc(
          list[str]: List of labels of SamplingFixture being managed.
          )pbdoc")
      .def_readonly("engine", &run_manager_type::engine,
                    R"pbdoc(
          libcasm.monte.RandomNumberEngine: Random number engine
          )pbdoc")
      .def_readonly("global_cutoff", &run_manager_type::engine,
                    R"pbdoc(
          bool: Global run cutoff if any sampling fixture is complete

          If true, the run is complete if any sampling fixture is complete.
          Otherwise, all sampling fixtures must be completed for the run to be
          completed.
          )pbdoc")
      .def(
          "sampling_fixture",
          [](run_manager_type const &self,
             std::string label) -> std::shared_ptr<sampling_fixture_type> {
            for (auto &fixture_ptr : self.sampling_fixtures) {
              if (fixture_ptr->label() == label) {
                return fixture_ptr;
              }
            }
            std::stringstream msg;
            msg << "Error in RunManager.sampling_fixture: label=" << label
                << " does not exist.";
            throw std::runtime_error(msg.str());
          },
          R"pbdoc(
          Get sampling fixture by label

          Parameters
          ----------
          label: str
              Label of sampling fixture to return.

          Returns
          -------
          sampling_fixture: SamplingFixture
              Sampling fixture with matching label.
          )pbdoc",
          py::arg("label"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
