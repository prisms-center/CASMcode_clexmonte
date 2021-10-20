#include "casm/clexmonte/canonical/sampling_functions.hh"

#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/sampling_functions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param system_data Reference to system_type, which can be used by
///     sampling functions to access the cluster expansion and composition
///     axes. The lifetime of the object it references must be
///     greater than the sampling functions that use the reference.
///
monte::StateSamplingFunctionMap<Configuration> make_sampling_functions(
    std::shared_ptr<system_type> const &system_data) {
  std::vector<monte::StateSamplingFunction<Configuration>> functions = {
      make_temperature_f(system_data), make_comp_n_f(system_data),
      make_comp_x_f(system_data), make_formation_energy_corr_f(system_data),
      make_formation_energy_f(system_data)};

  monte::StateSamplingFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
};

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
