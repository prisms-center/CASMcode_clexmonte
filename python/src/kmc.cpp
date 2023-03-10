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

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_kmc, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Kinetic Monte Carlo classes and methods

        libcasm.kmc
        -----------

        The libcasm.kmc package contains data structures and methods for running kinetic Monte Carlo calculations.

    )pbdoc";

  m.def(
      "placeholder", [](int a, int b) { return a + b; },
      R"pbdoc(
      Placeholder for docs generation

      Parameters
      ----------
      a: int
          LHS integer.

      b: int
          RHS integer.

      Returns
      -------
      sum: int
          The sum of a and b.
        )pbdoc",
      py::arg("a"), py::arg("b"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
