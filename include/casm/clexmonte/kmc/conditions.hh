#ifndef CASM_clexmonte_kmc_conditions
#define CASM_clexmonte_kmc_conditions

#include "casm/monte/state/ValueMap.hh"

namespace CASM {
namespace clexmonte {
namespace kmc {

/// \brief Holds conditions to avoid map lookup, including beta
struct Conditions {
  explicit Conditions(monte::ValueMap _conditions)
      : conditions(_conditions),
        beta(1.0 / (CASM::KB * conditions.scalar_values.at("temperature"))) {}

  monte::ValueMap conditions;
  double beta;
};

}  // namespace kmc
}  // namespace clexmonte
}  // namespace CASM

#endif
