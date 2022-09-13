#ifndef CASM_clexmonte_definitions
#define CASM_clexmonte_definitions

#include <map>
#include <string>

#include "casm/monte/definitions.hh"

namespace CASM {

namespace composition {
class CompositionCalculator;
class CompositionConverter;
}  // namespace composition

namespace clexmonte {
struct Configuration;
struct System;
}  // namespace clexmonte

namespace clexmonte {

typedef System system_type;

typedef Configuration config_type;

typedef monte::State<config_type> state_type;

// generic ConfigGenerator, and supported implementations
typedef monte::ConfigGenerator<config_type, state_type> config_generator_type;

// generic StateGenerator, and supported implementations
typedef monte::StateGenerator<config_type, state_type> state_generator_type;

typedef monte::Results<config_type> results_type;

// generic ResultsIO, and supported implementations
typedef monte::ResultsIO<config_type> results_io_type;

typedef std::function<Eigen::VectorXd(std::vector<Eigen::VectorXd> const &)>
    CorrCalculatorFunction;

}  // namespace clexmonte
}  // namespace CASM

#endif
