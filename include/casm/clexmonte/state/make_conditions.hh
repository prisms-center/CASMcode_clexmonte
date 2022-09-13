#ifndef CASM_clexmonte_state_make_conditions
#define CASM_clexmonte_state_make_conditions

#include <map>
#include <string>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/eigen.hh"

namespace CASM {

namespace composition {
class CompositionConverter;
}

namespace clexmonte {

// --- Composition ---

/// \brief Helper for making a conditions ValueMap, mol_composition
Eigen::VectorXd make_mol_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp, bool do_not_convert = false);

/// \brief Helper for making a conditions ValueMap, mol_composition
/// increment
Eigen::VectorXd make_mol_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp, bool do_not_convert = false);

/// \brief Helper for making a conditions ValueMap, param_composition
Eigen::VectorXd make_param_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp, bool do_not_convert = false);

/// \brief Helper for making a conditions ValueMap, param_composition
/// increment
Eigen::VectorXd make_param_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp, bool do_not_convert = false);

// --- Chemical potential ---

/// \brief Helper for making a conditions ValueMap, param_chem_pot
Eigen::VectorXd make_param_chem_pot(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

/// \brief Helper for making a conditions ValueMap, param_chem_pot
///     increment
Eigen::VectorXd make_param_chem_pot_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

}  // namespace clexmonte
}  // namespace CASM

#endif
