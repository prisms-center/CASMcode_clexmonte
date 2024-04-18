
#include "casm/clexmonte/monte_calculator/BaseMonteCalculator.hh"

namespace CASM {
namespace clexmonte {

BaseMonteCalculator::BaseMonteCalculator(
    std::string _calculator_name, std::set<std::string> _required_basis_set,
    std::set<std::string> _required_local_basis_set,
    std::set<std::string> _required_clex,
    std::set<std::string> _required_multiclex,
    std::set<std::string> _required_local_clex,
    std::set<std::string> _required_local_multiclex,
    std::set<std::string> _required_dof_spaces,
    std::set<std::string> _required_params,
    std::set<std::string> _optional_params, bool _time_sampling_allowed,
    bool _update_species, bool _is_multistate_method)
    : calculator_name(_calculator_name),
      required_basis_set(_required_basis_set),
      required_local_basis_set(_required_local_basis_set),
      required_clex(_required_clex),
      required_multiclex(_required_multiclex),
      required_local_clex(_required_local_clex),
      required_local_multiclex(_required_local_multiclex),
      required_dof_spaces(_required_dof_spaces),
      required_params(_required_params),
      optional_params(_optional_params),
      time_sampling_allowed(_time_sampling_allowed),
      update_species(_update_species),
      is_multistate_method(_is_multistate_method) {}

BaseMonteCalculator::~BaseMonteCalculator() {}

/// \brief Standardized check for whether system has required data
void BaseMonteCalculator::_check_system() const {
  for (std::string const &key : this->required_basis_set) {
    if (!is_basis_set(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' basis_set.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_basis_set) {
    if (!is_local_basis_set(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_basis_set.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_clex) {
    if (!is_clex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' clex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_multiclex) {
    if (!is_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' multiclex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_clex) {
    if (!is_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_clex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_local_multiclex) {
    if (!is_local_multiclex_data(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' local_multiclex.";
      throw std::runtime_error(msg.str());
    }
  }
  for (std::string const &key : this->required_dof_spaces) {
    if (!is_dof_space(*this->system, key)) {
      std::stringstream msg;
      msg << "Error preparing MonteCalculator " << this->calculator_name
          << ": no '" << key << "' dof_space.";
      throw std::runtime_error(msg.str());
    }
  }
}

/// \brief Validate calculator_params (for top-level key existence only)
///
/// Notes:
/// - Keys that start with "_" are ignored
/// - A warning is printed for keys that are not in the
///   required_params or optional_params sets
///
void BaseMonteCalculator::_check_params() const {
  for (auto key : this->required_params) {
    if (!params.contains(key)) {
      std::stringstream msg;
      msg << "Error: Missing required parameter '" << key << "'.";
      throw std::runtime_error(msg.str());
    }
  }

  auto &log = CASM::log();
  for (auto it = params.begin(); it != params.end(); ++it) {
    std::string key = it.name();
    if (key.empty()) {
      std::stringstream msg;
      msg << "Error: Empty parameter key.";
      throw std::runtime_error(msg.str());
    }
    if (key[0] == '_') {
      continue;
    }
    if (!this->required_params.count(key) &&
        !this->optional_params.count(key)) {
      log.indent() << "Warning: Unknown parameter '" << key << "'."
                   << std::endl;
    }
  }
}

/// \brief Clone the BaseMonteCalculator
std::unique_ptr<BaseMonteCalculator> BaseMonteCalculator::clone() const {
  return std::unique_ptr<BaseMonteCalculator>(this->_clone());
}

}  // namespace clexmonte
}  // namespace CASM
