#ifndef CASM_clexmonte_canonical_definitions
#define CASM_clexmonte_canonical_definitions

#include <map>
#include <string>

class MTRand;

namespace CASM {

namespace monte {

struct CompletionCheckParams;
struct SamplingParams;

template <typename _ConfigType, typename _RunInfoType>
class ConfigGenerator;
template <typename _ConfigType>
class FixedConfigGenerator;

template <typename _ConfigType>
struct Results;

template <typename _ConfigType>
class ResultsIO;

template <typename _ConfigType>
struct State;

template <typename _ConfigType, typename _RunInfoType>
class StateGenerator;
template <typename _ConfigType>
class IncrementalConditionsStateGenerator;

template <typename _ConfigType>
class StateSampler;

template <typename _ConfigType>
struct StateSamplingFunction;

template <typename _ConfigType>
using StateSamplingFunctionMap =
    std::map<std::string, StateSamplingFunction<_ConfigType>>;

}  // namespace monte

namespace clexmonte {

struct Configuration;

struct OccSystem;

}  // namespace clexmonte

namespace clexmonte {
namespace canonical {

struct canonical_tag {};

typedef OccSystem system_type;

typedef Configuration config_type;

typedef monte::State<config_type> state_type;

// generic ConfigGenerator, and supported implementations
typedef monte::ConfigGenerator<config_type, state_type> config_generator_type;
typedef monte::FixedConfigGenerator<config_type> fixed_config_generator_type;

// generic StateGenerator, and supported implementations
typedef monte::StateGenerator<config_type, state_type> state_generator_type;
typedef monte::IncrementalConditionsStateGenerator<config_type>
    incremental_state_generator_type;

typedef std::map<std::string, monte::StateSampler<config_type>>
    samplers_map_type;

typedef monte::Results<config_type> results_type;
typedef monte::ResultsIO<config_type> results_io_type;

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
