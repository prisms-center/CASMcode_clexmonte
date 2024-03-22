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
#include "casm/clexmonte/semigrand_canonical/calculator_impl.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexmonte/semigrand_canonical/json_io.hh"
#include "casm/clexmonte/semigrand_canonical/potential.hh"
#include "casm/monte/RandomNumberGenerator.hh"

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
typedef clexmonte::System system_type;
typedef clexmonte::run_manager_type<engine_type> run_manager_type;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_semigrand_canonical, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Cluster expansion semi-grand canonical Monte Carlo

        libcasm.clexmonte.semigrand_canonical._clexmonte_semigrand_canonical
        --------------------------------------------------------------------

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
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.events");
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
                    &conditions_type::temperature_is_zero_tol)
      .def_readonly("temperature", &conditions_type::temperature)
      .def_readonly("beta", &conditions_type::beta)
      .def("set_temperature",
           [](conditions_type &self, double temperature) {
             self.set_temperature(temperature);
           })
      .def("set_temperature_from_value_map",
           [](conditions_type &self, monte::ValueMap const &map) {
             self.set_temperature(map);
           })
      .def("put_temperature", &conditions_type::put_temperature)
      .def_readonly("composition_converter",
                    &conditions_type::composition_converter)
      .def_readonly("param_chem_pot", &conditions_type::param_chem_pot)
      .def_readonly("exchange_chem_pot", &conditions_type::exchange_chem_pot)
      .def("set_param_chem_pot",
           [](conditions_type &self, Eigen::VectorXd const &param_chem_pot) {
             self.set_param_chem_pot(param_chem_pot);
           })
      .def("set_param_chem_pot_from_value_map",
           [](conditions_type &self, monte::ValueMap const &map) {
             self.set_param_chem_pot(map);
           })
      .def("put_param_chem_pot", &conditions_type::put_param_chem_pot)
      .def("set_all", &conditions_type::set_all)
      .def("to_value_map", &conditions_type::to_value_map);

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
      .def("set", &potential_type::set)
      .def("state", &potential_type::state)
      .def("conditions", &potential_type::conditions)
      .def("formation_energy", &potential_type::formation_energy)
      .def("per_supercell", &potential_type::per_supercell)
      .def("per_unitcell", &potential_type::per_unitcell)
      .def("occ_delta_per_supercell", &potential_type::occ_delta_per_supercell);

  py::class_<calculator_type, std::shared_ptr<calculator_type>>(
      m, "SemiGrandCanonicalCalculator",
      R"pbdoc(
      Implements semi-grand canonical Monte Carlo calculations
      )pbdoc")
      .def(py::init<std::shared_ptr<system_type>>(),
           R"pbdoc(
        .. rubric:: Constructor

        Parameters
        ----------
        system : libcasm.clexmonte.System
            Cluster expansion model system data.
        )pbdoc",
           py::arg("system"))
      .def(
          "run",  //&calculator_type::run,
          [](calculator_type &self, state_type &state,
             monte::OccLocation *occ_location, run_manager_type &run_manager) {
            if (!occ_location) {
              monte::Conversions const &convert =
                  get_index_conversions(*self.system, state);
              monte::OccCandidateList const &occ_candidate_list =
                  get_occ_candidate_list(*self.system, state);
              monte::OccLocation tmp(convert, occ_candidate_list,
                                     self.update_species);
              tmp.initialize(get_occupation(state));
              self.run(state, tmp, run_manager);
            } else {
              self.run(state, *occ_location, run_manager);
            }
          },
          R"pbdoc(
          Perform a single run, evolving the input state

          Parameters
          ----------
          state : libcasm.clexmonte.State
              The input state.
          run_manager: libcasm.clexmonte.RunManager
              Specifies sampling and convergence criteria and collects results
          occ_location: Optional[libcasm.monte.events.OccLocation] = None
              Current occupant location list. If provided, the user is
              responsible for ensuring it is up-to-date with the current
              occupation of `state` and it is used and updated during the run.
              If None, a occupant location list is generated for the run.
          )pbdoc",
          py::arg("state"), py::arg("run_manager"), py::arg("occ_location"))
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
