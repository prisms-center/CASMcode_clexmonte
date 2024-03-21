#include "KMCCompleteEventListTestSystem.hh"

#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/events/OccLocation.hh"

namespace test {
KMCCompleteEventListTestSystem::KMCCompleteEventListTestSystem()
    : KMCTestSystem() {}

KMCCompleteEventListTestSystem::KMCCompleteEventListTestSystem(
    std::string _project_name, std::string _test_dir_name,
    fs::path _input_file_path)
    : KMCTestSystem(_project_name, _test_dir_name, _input_file_path) {}

void KMCCompleteEventListTestSystem::make_prim_event_list(
    std::vector<std::string> const &clex_names,
    std::vector<std::string> const &multiclex_names) {
  prim_event_list = clexmonte::make_prim_event_list(*system);
  prim_impact_info_list = clexmonte::make_prim_impact_info_list(
      *system, prim_event_list, clex_names, multiclex_names);
}

/// \brief Make complete event list
///
/// Note: This calls occ_location->initialize. For correct atom tracking and
/// stochastic canonical / semi-grand canoncical event choosing,
/// occ_location->initialize must be called again if the configuration is
/// modified directly instead of via occ_location->apply. Event calculations
/// would be still be correct.
void KMCCompleteEventListTestSystem::make_complete_event_list(
    monte::State<clexmonte::Configuration> const &state) {
  occ_location = std::make_unique<monte::OccLocation>(
      get_index_conversions(*system, state),
      get_occ_candidate_list(*system, state));
  occ_location->initialize(get_occupation(state));

  event_list = clexmonte::make_complete_event_list(
      prim_event_list, prim_impact_info_list, *occ_location);
}

}  // namespace test
