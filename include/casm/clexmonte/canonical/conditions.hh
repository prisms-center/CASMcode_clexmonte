#ifndef CASM_clexmonte_canonical_conditions
#define CASM_clexmonte_canonical_conditions

#include "casm/clexmonte/canonical/definitions.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

/// \brief Helper for making a conditions VectorValueMap for canonical Monte
///     Carlo calculations
monte::VectorValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions VectorValueMap for canonical Monte
///     Carlo calculations
monte::VectorValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
