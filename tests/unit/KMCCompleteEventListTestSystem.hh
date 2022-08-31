#ifndef CASM_unittest_KMCCompleteEventListTestSystem
#define CASM_unittest_KMCCompleteEventListTestSystem

#include "KMCTestSystem.hh"
#include "casm/clexmonte/events/event_methods.hh"
#include "casm/clexmonte/kmc/CompleteEventList.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/events/OccLocation.hh"

using namespace CASM;

class KMCCompleteEventListTestSystem : public KMCTestSystem {
 public:
  std::vector<clexmonte::PrimEventData> prim_event_list;
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  std::unique_ptr<monte::OccLocation> occ_location;
  clexmonte::kmc::CompleteEventList event_list;

  KMCCompleteEventListTestSystem() : KMCTestSystem() {}

  KMCCompleteEventListTestSystem(std::string _project_name,
                                 std::string _test_dir_name,
                                 fs::path _input_file_path)
      : KMCTestSystem(_project_name, _test_dir_name, _input_file_path) {}

  void make_prim_event_list(
      std::vector<std::string> const &clex_names = {"formation_energy"},
      std::vector<std::string> const &multiclex_names = {}) {
    prim_event_list = clexmonte::make_prim_event_list(*system);
    prim_impact_info_list = clexmonte::make_prim_impact_info_list(
        *system, prim_event_list, clex_names, multiclex_names);
  }

  void make_complete_event_list(
      monte::State<clexmonte::Configuration> const &state) {
    occ_location = std::make_unique<monte::OccLocation>(
        get_index_conversions(*system, state),
        get_occ_candidate_list(*system, state));
    occ_location->initialize(get_occupation(state));

    event_list = clexmonte::kmc::make_complete_event_list(
        prim_event_list, prim_impact_info_list, *occ_location);
  }
};

#endif
