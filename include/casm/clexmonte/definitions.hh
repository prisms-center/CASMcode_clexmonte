#ifndef CASM_clexmonte_definitions
#define CASM_clexmonte_definitions

#include <map>
#include <string>

#include "casm/monte/BasicStatistics.hh"
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

typedef monte::BasicStatistics statistics_type;

typedef monte::State<config_type> state_type;

typedef monte::StateSamplingFunction<config_type> state_sampling_function_type;
typedef monte::ResultsAnalysisFunction<config_type, statistics_type>
    results_analysis_function_type;
typedef monte::StateModifyingFunction<config_type>
    state_modifying_function_type;
typedef monte::RunManager<config_type, statistics_type> run_manager_type;
typedef monte::SamplingFixtureParams<config_type, statistics_type>
    sampling_figure_params_type;

struct Conditions;

// generic ConfigGenerator, and supported implementations
typedef monte::ConfigGenerator<config_type, monte::RunData<config_type>>
    config_generator_type;

// generic StateGenerator, and supported implementations
typedef monte::StateGenerator<config_type, monte::RunData<config_type>>
    state_generator_type;

typedef monte::RunData<config_type> run_data_type;
typedef monte::Results<config_type, statistics_type> results_type;

// generic ResultsIO, and supported implementations
typedef monte::ResultsIO<results_type> results_io_type;

typedef std::function<Eigen::VectorXd(std::vector<Eigen::VectorXd> const &)>
    CorrCalculatorFunction;

}  // namespace clexmonte
}  // namespace CASM

#endif
