#ifndef CASM_clexmonte_system_parse_composition
#define CASM_clexmonte_system_parse_composition

#include "casm/global/eigen.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace composition {
class CompositionConverter;
}

namespace clexmonte {

/// \brief Parse composition (either comp_n or comp_x) from JSON
void parse_composition_as_provided(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool &is_comp_n);

/// \brief Parse composition from JSON, return as comp_n or comp_n increment
void parse_composition_for_comp_n(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment);

/// \brief Parse composition from JSON, return as comp_x or comp_x increment
void parse_composition_for_comp_x(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment);

}  // namespace clexmonte
}  // namespace CASM

#endif
