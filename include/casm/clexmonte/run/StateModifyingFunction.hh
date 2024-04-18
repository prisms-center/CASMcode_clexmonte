#ifndef CASM_clexmonte_StateModifyingFunction
#define CASM_clexmonte_StateModifyingFunction

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/definitions.hh"

namespace CASM {
namespace clexmonte {

struct StateModifyingFunction {
  /// \brief Constructor - default component names
  StateModifyingFunction(
      std::string _name, std::string _description,
      std::function<void(state_type &, monte::OccLocation *)> _function)
      : name(_name), description(_description), function(_function) {}

  std::string name;

  std::string description;

  std::function<void(state_type &, monte::OccLocation *)> function;

  /// \brief Evaluates `function`
  void operator()(state_type &state, monte::OccLocation *occ_location) const {
    function(state, occ_location);
  }
};

}  // namespace clexmonte
}  // namespace CASM

#endif
