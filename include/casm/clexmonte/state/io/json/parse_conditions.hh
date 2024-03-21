#ifndef CASM_clexmonte_state_parse_conditions
#define CASM_clexmonte_state_parse_conditions

#include <memory>

#include "casm/clexmonte/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

void parse(InputParser<Conditions> &parser,
           std::shared_ptr<system_type> const &system, bool is_increment);

/// \brief Parse all conditions
void parse_conditions(InputParser<monte::ValueMap> &parser,
                      std::shared_ptr<system_type> const &system,
                      bool is_increment);

/// \brief Parse temperature scalar value
void parse_temperature(InputParser<monte::ValueMap> &parser);

/// \brief Parse "mol_composition" and store as
///     "mol_composition" vector values
void parse_mol_composition(InputParser<monte::ValueMap> &parser,
                           std::shared_ptr<system_type> const &system,
                           bool is_increment);

/// \brief Parse "param_composition" and store as
///     "param_composition" vector values
void parse_param_composition(InputParser<monte::ValueMap> &parser,
                             std::shared_ptr<system_type> const &system,
                             bool is_increment);

/// \brief Parse "param_chem_pot" and store as
///     "param_chem_pot" vector values
void parse_param_chem_pot(InputParser<monte::ValueMap> &parser,
                          std::shared_ptr<system_type> const &system);

/// \brief Parse "include_formation_energy" or set default value
void parse_include_formation_energy(InputParser<monte::ValueMap> &parser);

/// \brief Parse boolean value to monte::ValueMap
void parse_boolean(InputParser<monte::ValueMap> &parser, std::string option);

/// \brief Parse scalar value to monte::ValueMap
void parse_scalar(InputParser<monte::ValueMap> &parser, std::string option);

/// \brief Parse vector value to monte::ValueMap
void parse_vector(InputParser<monte::ValueMap> &parser, std::string option);

/// \brief Parse matrix value to monte::ValueMap
void parse_matrix(InputParser<monte::ValueMap> &parser, std::string option);

/// \brief Parse "corr_matching_pot"
void parse_corr_matching_pot(InputParser<monte::ValueMap> &parser,
                             bool is_increment);

/// \brief Parse "random_alloy_corr_matching_pot"
void parse_random_alloy_corr_matching_pot(
    InputParser<monte::ValueMap> &parser,
    std::shared_ptr<system_type> const &system, bool is_increment);

// ~~~ Templated, ConditionsType parsing ~~~

/// \brief Parse 'temperature' (as a required attribute)
template <typename ConditionsType>
void parse_temperature(InputParser<ConditionsType> &parser);

/// \brief Parse 'param_chem_pot' (as a required attribute)
template <typename ConditionsType>
void parse_param_chem_pot(
    InputParser<ConditionsType> &parser,
    composition::CompositionConverter const &composition_converter);

/// \brief Parse 'mol_composition' and/or 'param_composition' (at least one of
/// which is required, if both are provided they must be equivalent)
template <typename ConditionsType>
void parse_composition(
    InputParser<ConditionsType> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment);
}  // namespace clexmonte
}  // namespace CASM

#endif
