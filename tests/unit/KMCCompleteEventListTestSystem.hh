#ifndef CASM_unittest_KMCCompleteEventListTestSystem
#define CASM_unittest_KMCCompleteEventListTestSystem

#include "KMCTestSystem.hh"
#include "casm/clexmonte/events/CompleteEventList.hh"
#include "casm/monte/events/OccLocation.hh"

namespace test {

using namespace CASM;

class KMCCompleteEventListTestSystem : public KMCTestSystem {
 public:
  std::vector<clexmonte::PrimEventData> prim_event_list;
  std::vector<clexmonte::EventImpactInfo> prim_impact_info_list;

  std::unique_ptr<monte::OccLocation> occ_location;
  clexmonte::CompleteEventList event_list;

  KMCCompleteEventListTestSystem();

  KMCCompleteEventListTestSystem(std::string _project_name,
                                 std::string _test_dir_name,
                                 fs::path _input_file_path);

  void make_prim_event_list(
      std::vector<std::string> const &clex_names = {"formation_energy"},
      std::vector<std::string> const &multiclex_names = {});

  /// \brief Make complete event list
  ///
  /// Note: This calls occ_location->initialize. For correct atom tracking and
  /// stochastic canonical / semi-grand canoncical event choosing,
  /// occ_location->initialize must be called again if the configuration is
  /// modified directly instead of via occ_location->apply. Event calculations
  /// would be still be correct.
  void make_complete_event_list(
      monte::State<clexmonte::Configuration> const &state);
};

}  // namespace test
#endif
