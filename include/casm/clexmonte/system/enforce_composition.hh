#ifndef CASM_clexmonte_system_enforce_composition
#define CASM_clexmonte_system_enforce_composition

#include <vector>

#include "casm/global/eigen.hh"

class MTRand;

namespace CASM {

namespace composition {
class CompositionCalculator;
}

namespace monte {
class Conversions;
class OccSwap;
}  // namespace monte

namespace clexmonte {

/// \brief Apply grand canonical swaps to enforce target composition
void enforce_composition(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_comp_n,
    composition::CompositionCalculator const &composition_calculator,
    monte::Conversions const &convert,
    std::vector<monte::OccSwap> const &grand_canonical_swaps,
    MTRand &random_number_generator);

}  // namespace clexmonte
}  // namespace CASM

#endif
