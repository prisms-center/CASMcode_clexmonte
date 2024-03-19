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
#include "casm/clexmonte/semigrand_canonical"
#include "casm/monte/RandomNumberGenerator.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// used for libcasm.clexmonte:
typedef std::mt19937_64 engine_type;
typedef RandomNumberGenerator<engine_type> generator_type;
// clexmonte::system_type
// clexmonte::config_type
// clexmonte::state_type

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_semigrand_canonical, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Cluster expansion semi-grand canonical Monte Carlo

        libcasm.clexmonte._semigrand_canonical
        --------------------------------------

        Includes:

        - the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalConditions`
          class for representing thermodynamic conditions,
        - the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalEventGenerator`
          class for proposing events in the semi-grand canonical ensemble,
        - the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalPotential`
          class for calculating changes in the semi-grand canonical energy due to the
          proposed events, and
        - the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalCalculator`
          class for sampling microstates in the semi-grand canonical ensemble.

        It also makes use of:

        - :class:`~libcasm.monte.sampling.SamplingFixture` and
          :class:`~libcasm.monte.sampling.RunManager`, to control sampling, convergence
          checking, and results output.


    )pbdoc";
  //  py::module::import("libcasm.clexulator");
  //  py::module::import("libcasm.composition");
  //  py::module::import("libcasm.configuration");
  //  py::module::import("libcasm.xtal");
  py::module::import("libcasm.monte");
  py::module::import("libcasm.clexmonte._system");
  py::module::import("libcasm.clexmonte._state");

  py::class_<
      clexmonte::semigrand_canonical::SemiGrandCanonical,
      std::shared_ptr<clexmonte::semigrand_canonical::SemiGrandCanonical>>(
      m, "SemiGrandCanonicalCalculator",
      R"pbdoc(
      Implements semi-grand canonical Monte Carlo calculations
      )pbdoc")
      .def(py::init<std::shared_ptr<system_type>,
                    std::shared_ptr<engine_type>>(),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        system : libcasm.clexmonte.System
            Cluster expansion model system data.
        engine : Optional[:class:`~libcasm.monte.RandomNumberEngine`] = None
              A :class:`~libcasm.monte.RandomNumberEngine` to use for
              generating random numbers. If provided, the engine will be
              shared. If None, then a new
              :class:`~libcasm.monte.RandomNumberEngine` will be constructed
              and seeded using std::random_device.
        )pbdoc",
           py::arg("system"), py::arg("engine") = py::none());

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
