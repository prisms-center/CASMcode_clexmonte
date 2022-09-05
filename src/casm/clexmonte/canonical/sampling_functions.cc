#include "casm/clexmonte/canonical/sampling_functions.hh"

#include "casm/clexmonte/state/sampling_functions.hh"
#include "casm/clexmonte/system/System.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Construct functions that may be used to sample various quantities of
///     the Monte Carlo calculation as it runs
///
/// \param system Shared pointer to System data, which can be used by
///     sampling functions to access data such as the prim, the cluster
///     expansion, and the composition axes.
/// \param tag The canonical_tag is used to help overload disambiguation
///
monte::StateSamplingFunctionMap<Configuration> make_sampling_functions(
    std::shared_ptr<system_type> const &system, canonical_tag tag) {
  std::vector<monte::StateSamplingFunction<Configuration>> functions = {
      make_temperature_f(system),       make_mol_composition_f(system),
      make_param_composition_f(system), make_formation_energy_corr_f(system),
      make_formation_energy_f(system),  make_potential_energy_f(system)};

  monte::StateSamplingFunctionMap<Configuration> function_map;
  for (auto const &f : functions) {
    function_map.emplace(f.name, f);
  }
  return function_map;
};

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM
