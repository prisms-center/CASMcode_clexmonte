#ifndef CASM_clexmonte_StateModifyingFunction
#define CASM_clexmonte_StateModifyingFunction

#include "casm/monte/definitions.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

struct StateModifyingFunction {
  /// \brief Constructor - default component names
  StateModifyingFunction(std::string _name, std::string _description,
                         std::function<void(state_type &)> _function)
      : name(_name), description(_description), function(_function) {}

  std::string name;

  std::string description;

  std::function<void(state_type &)> function;

  /// \brief Evaluates `function`
  void operator()(state_type &state) const { function(state); }
};

}  // namespace clexmonte
}  // namespace CASM

#endif
