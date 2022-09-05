#ifndef CASM_clexmonte_state_enforce_composition
#define CASM_clexmonte_state_enforce_composition

#include <vector>

#include "casm/global/eigen.hh"

class MTRand;

namespace CASM {

namespace composition {
class CompositionCalculator;
}

namespace monte {
class OccSwap;
class OccLocation;
}  // namespace monte

namespace clexmonte {

/// \brief Apply grand canonical swaps to enforce target composition
void enforce_composition(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_mol_composition,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<monte::OccSwap> const &grand_canonical_swaps,
    monte::OccLocation &occ_location, MTRand &random_number_generator);

}  // namespace clexmonte
}  // namespace CASM

#endif
