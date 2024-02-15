#ifndef CASM_unittest_KMCCompleteEventCalculatorTestSystem
#define CASM_unittest_KMCCompleteEventCalculatorTestSystem

#include "KMCCompleteEventListTestSystem.hh"
#include "casm/clexmonte/kinetic/kinetic_events.hh"

namespace test {

using namespace CASM;

class KMCCompleteEventCalculatorTestSystem
    : public KMCCompleteEventListTestSystem {
 public:
  std::shared_ptr<clexmonte::Conditions> conditions;
  std::vector<clexmonte::kinetic::EventStateCalculator> prim_event_calculators;
  std::shared_ptr<clexmonte::kinetic::CompleteEventCalculator> event_calculator;

  KMCCompleteEventCalculatorTestSystem();

  KMCCompleteEventCalculatorTestSystem(std::string _project_name,
                                       std::string _test_dir_name,
                                       fs::path _input_file_path);

  // Note: For correct atom tracking and stochastic canonical / semi-grand
  // canoncical
  //  event choosing, after this, occ_location->initialize must be called again
  // if the configuration is modified directly instead of via
  // occ_location->apply. Event calculations would be still be correct.
  void make_complete_event_calculator(
      monte::State<clexmonte::Configuration> const &state,
      std::vector<std::string> const &clex_names = {"formation_energy"},
      std::vector<std::string> const &multiclex_names = {});
};
}  // namespace test
#endif
