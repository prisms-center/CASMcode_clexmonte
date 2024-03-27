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
#include "casm/clexmonte/state/enforce_composition.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccLocation.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;
typedef std::mt19937_64 engine_type;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clexmonte_functions, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
      Cluster expansion Monte Carlo functions
      )pbdoc";
  py::module::import("libcasm.composition");
  py::module::import("libcasm.monte");
  py::module::import("libcasm.monte.events");

  m.def(
      "enforce_composition",
      [](clexmonte::state_type &state,
         Eigen::VectorXd const &target_mol_composition,
         std::optional<std::shared_ptr<clexmonte::System>> _system,
         std::optional<composition::CompositionCalculator>
             composition_calculator,
         std::optional<std::vector<monte::OccSwap>> semigrand_canonical_swaps,
         monte::OccLocation *occ_location,
         std::optional<std::shared_ptr<engine_type>> _engine) {
        std::shared_ptr<clexmonte::System> system;
        if (_system.has_value()) {
          system = _system.value();
        }

        std::shared_ptr<engine_type> engine;
        if (_engine.has_value()) {
          engine = _engine.value();
        }

        std::unique_ptr<monte::OccLocation> _occ_location;
        if (!occ_location) {
          if (!system) {
            throw std::runtime_error(
                "Error in enforce_composition: "
                "occ_location is None and system is None");
          }
          monte::Conversions const &convert =
              get_index_conversions(*system, state);
          monte::OccCandidateList const &occ_candidate_list =
              get_occ_candidate_list(*system, state);
          bool update_species = false;
          _occ_location = std::make_unique<monte::OccLocation>(
              convert, occ_candidate_list, update_species);
          _occ_location->initialize(clexmonte::get_occupation(state));
          occ_location = _occ_location.get();
        }

        auto pick_cc = [=]() -> composition::CompositionCalculator const & {
          if (composition_calculator.has_value()) {
            return *composition_calculator;
          } else {
            if (!system) {
              throw std::runtime_error(
                  "Error in enforce_composition: "
                  "composition_calculator is None and system is None");
            }
            return system->composition_calculator;
          }
        };

        composition::CompositionCalculator const &cc = pick_cc();
        if (target_mol_composition.size() != cc.components().size()) {
          std::stringstream msg;
          msg << "Error in enforce_composition: mismatch ";
          msg << "target_mol_composition=["
              << target_mol_composition.transpose();
          msg << "], expected: [" << cc.components() << "]" << std::endl;
          throw std::runtime_error(msg.str());
        }

        auto pick_swaps = [=]() -> std::vector<monte::OccSwap> const & {
          if (semigrand_canonical_swaps.has_value()) {
            return *semigrand_canonical_swaps;
          } else {
            if (!system) {
              throw std::runtime_error(
                  "Error in enforce_composition: "
                  "semigrand_canonical_swaps is None and system is None");
            }
            return system->semigrand_canonical_swaps;
          }
        };

        monte::RandomNumberGenerator<engine_type> random_number_generator(
            engine);

        clexmonte::enforce_composition(
            get_occupation(state), target_mol_composition, pick_cc(),
            pick_swaps(), *occ_location, random_number_generator);
      },
      R"pbdoc(
            Apply grand canonical swaps to enforce a target composition

            .. rubric:: Method

            - Find which of the provided grand canonical swap types transforms
              the composition most closely to the target composition
            - If no swap can improve the composition, return
            - Propose and apply an event consistent with the found swap type
            - Repeat

            Parameters
            ----------
            state: MonteCarloState
                The state to modify to enforce a target composition
            target_mol_composition: np.ndarray
                The target mol composition per unit cell, :math:`\vec{n}`.
            system: Optional[System] = None
                System data. Used to get the composition calculator, get the
                allowed semi-grand canonical swaps, and construct occupant
                location list, unless they are provided explicitly.
            composition_calculator: Optional[libcasm.composition.CompositionCalculator] = None
                Composition calculator. If not provided, the `system`
                composition calculator is used. Raises if neither `system` nor
                `composition_calculator` are provided.
            semigrand_canonical_swaps: Optional[list[libcasm.monte.events.OccSwap]] = None
                Swaps to use to enforce composition. If not provided, the
                `system` semi-grand canonical swaps are used. Raises if
                neither `system` nor `composition_calculator` are provided.
            occ_location: Optional[libcasm.monte.events.OccLocation] = None
                Current occupant location list. If provided, the user is
                responsible for ensuring it is up-to-date with the current
                occupation of `state` and it is used and updated during the
                run. If None, a occupant location list is generated for the
                function. Raises if neither `system` nor `occ_location` are
                provided.
            engine: Optional[libcasm.monte.RandomNumberEngine] = None
                Optional random number engine to use.
            )pbdoc",
      py::arg("state"), py::arg("target_mol_composition"),
      py::arg("system") = std::nullopt,
      py::arg("composition_calculator") = std::nullopt,
      py::arg("semigrand_canonical_swaps") = std::nullopt,
      py::arg("occ_location") = static_cast<monte::OccLocation *>(nullptr),
      py::arg("engine") = std::nullopt);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
