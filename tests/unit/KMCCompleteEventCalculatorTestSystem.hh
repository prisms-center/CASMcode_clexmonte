#ifndef CASM_unittest_KMCCompleteEventCalculatorTestSystem
#define CASM_unittest_KMCCompleteEventCalculatorTestSystem

#include "KMCCompleteEventListTestSystem.hh"
#include "casm/clexmonte/events/CompleteEventCalculator.hh"

using namespace CASM;

class KMCCompleteEventCalculatorTestSystem
    : public KMCCompleteEventListTestSystem {
 public:
  std::shared_ptr<clexmonte::Conditions> conditions;
  std::vector<clexmonte::EventStateCalculator> prim_event_calculators;
  std::shared_ptr<clexmonte::CompleteEventCalculator> event_calculator;

  KMCCompleteEventCalculatorTestSystem() : KMCCompleteEventListTestSystem() {}

  KMCCompleteEventCalculatorTestSystem(std::string _project_name,
                                       std::string _test_dir_name,
                                       fs::path _input_file_path)
      : KMCCompleteEventListTestSystem(_project_name, _test_dir_name,
                                       _input_file_path) {}

  // Note: For correct atom tracking and stochastic canonical / grand canoncical
  //  event choosing, after this, occ_location->initialize must be called again
  // if the configuration is modified directly instead of via
  // occ_location->apply. Event calculations would be still be correct.
  void make_complete_event_calculator(
      monte::State<clexmonte::Configuration> const &state,
      std::vector<std::string> const &clex_names = {"formation_energy"},
      std::vector<std::string> const &multiclex_names = {}) {
    this->make_prim_event_list(clex_names, multiclex_names);

    // Note: For correct atom tracking and stochastic canonical / grand
    // canoncical
    //  event choosing, after this, occ_location->initialize must be called
    //  again
    // if the configuration is modified directly instead of via
    // occ_location->apply. Event calculations would be still be correct.
    this->make_complete_event_list(state);

    /// Make std::shared_ptr<clexmonte::Conditions> object from state.conditions
    conditions = std::make_shared<clexmonte::Conditions>(
        clexmonte::make_conditions_from_value_map(
            state.conditions, *system->prim->basicstructure,
            system->composition_converter, clexmonte::CorrCalculatorFunction(),
            CASM::TOL));

    prim_event_calculators = clexmonte::make_prim_event_calculators(
        *system, state, prim_event_list, conditions);

    // Construct CompleteEventCalculator
    event_calculator = std::make_shared<clexmonte::CompleteEventCalculator>(
        prim_event_list, prim_event_calculators, event_list.events);
  }
};

#endif
