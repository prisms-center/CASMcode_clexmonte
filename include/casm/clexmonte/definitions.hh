#ifndef CASM_clexmonte_definitions
#define CASM_clexmonte_definitions

#include <map>
#include <random>
#include <string>

#include "casm/monte/BasicStatistics.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/definitions.hh"

namespace CASM {

namespace composition {
class CompositionCalculator;
class CompositionConverter;
}  // namespace composition

namespace config {
struct Configuration;
}

namespace clexmonte {
typedef config::Configuration Configuration;
struct System;
}  // namespace clexmonte

namespace clexmonte {

typedef System system_type;
typedef config::Configuration config_type;
typedef monte::BasicStatistics statistics_type;
typedef monte::State<config_type> state_type;
// typedef std::mt19937_64 default_engine_type;
typedef monte::default_engine_type default_engine_type;
struct Conditions;

// ### Sampling ###

typedef monte::StateSamplingFunction state_sampling_function_type;
typedef monte::jsonStateSamplingFunction json_state_sampling_function_type;
typedef monte::ResultsAnalysisFunction<config_type, statistics_type>
    results_analysis_function_type;
typedef monte::SamplingFixtureParams<config_type, statistics_type>
    sampling_fixture_params_type;

template <typename EngineType>
using run_manager_type =
    monte::RunManager<config_type, statistics_type, EngineType>;
typedef monte::Results<config_type, statistics_type> results_type;
typedef monte::ResultsIO<results_type> results_io_type;

// ### State generation ###

struct RunData;

class StateGenerator;
typedef StateGenerator state_generator_type;
class IncrementalConditionsStateGenerator;

class ConfigGenerator;
typedef ConfigGenerator config_generator_type;
class FixedConfigGenerator;

struct StateModifyingFunction;
using StateModifyingFunctionMap = std::map<std::string, StateModifyingFunction>;

typedef std::function<Eigen::VectorXd(std::vector<Eigen::VectorXd> const &)>
    CorrCalculatorFunction;

}  // namespace clexmonte
}  // namespace CASM

#endif
