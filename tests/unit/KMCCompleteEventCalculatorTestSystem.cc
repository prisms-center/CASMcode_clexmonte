
#include "KMCCompleteEventCalculatorTestSystem.hh"

namespace test {

KMCCompleteEventCalculatorTestSystem::KMCCompleteEventCalculatorTestSystem()
    : KMCCompleteEventListTestSystem() {}

KMCCompleteEventCalculatorTestSystem::KMCCompleteEventCalculatorTestSystem(
    std::string _project_name, std::string _test_dir_name,
    fs::path _input_file_path)
    : KMCCompleteEventListTestSystem(_project_name, _test_dir_name,
                                     _input_file_path) {}

// Note: For correct atom tracking and stochastic canonical / semi-grand
// canoncical
//  event choosing, after this, occ_location->initialize must be called again
// if the configuration is modified directly instead of via
// occ_location->apply. Event calculations would be still be correct.
void KMCCompleteEventCalculatorTestSystem::make_complete_event_calculator(
    monte::State<clexmonte::Configuration> const &state,
    std::vector<std::string> const &clex_names,
    std::vector<std::string> const &multiclex_names) {
  using namespace clexmonte;

  this->make_prim_event_list(clex_names, multiclex_names);

  // Note: For correct atom tracking and stochastic canonical / semi-grand
  // canoncical
  //  event choosing, after this, occ_location->initialize must be called
  //  again
  // if the configuration is modified directly instead of via
  // occ_location->apply. Event calculations would be still be correct.
  this->make_complete_event_list(state);

  /// Make std::shared_ptr<clexmonte::Conditions> object from state.conditions
  conditions = make_conditions(*system, state);

  prim_event_calculators = kinetic::make_prim_event_calculators(
      system, state, prim_event_list, conditions);

  // Construct CompleteEventCalculator
  event_calculator = std::make_shared<kinetic::CompleteEventCalculator>(
      prim_event_list, prim_event_calculators, event_list.events);
}

}  // namespace test
