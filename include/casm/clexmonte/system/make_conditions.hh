#ifndef CASM_clexmonte_system_make_conditions
#define CASM_clexmonte_system_make_conditions

#include <map>
#include <string>

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/monte/definitions.hh"

namespace CASM {

namespace composition {
class CompositionConverter;
}

namespace clexmonte {

// Note: For double (ex. "temperature") or vector<double> to Eigen::VectorXd
// - Use `Eigen::VectorXd to_VectorXd(double value)`
// - Use `Eigen::VectorXd to_VectorXd(std::vector<double> const &value)`

// --- Composition ---

/// \brief Helper for making a conditions VectorValueMap, mol_composition
Eigen::VectorXd make_mol_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions VectorValueMap, mol_composition
/// increment
Eigen::VectorXd make_mol_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions VectorValueMap, param_composition
Eigen::VectorXd make_param_composition(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

/// \brief Helper for making a conditions VectorValueMap, param_composition
/// increment
Eigen::VectorXd make_param_composition_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp);

// --- Chemical potential ---

/// \brief Helper for making a conditions VectorValueMap, mol_chem_pot
Eigen::VectorXd make_mol_chem_pot(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

/// \brief Helper for making a conditions VectorValueMap, mol_chem_pot increment
Eigen::VectorXd make_mol_chem_pot_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

/// \brief Helper for making a conditions VectorValueMap, param_chem_pot
Eigen::VectorXd make_param_chem_pot(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

/// \brief Helper for making a conditions VectorValueMap, param_chem_pot
///     increment
Eigen::VectorXd make_param_chem_pot_increment(
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> chem_pot);

}  // namespace clexmonte
}  // namespace CASM

#endif
