#ifndef CASM_clexmonte_run_covariance_functions
#define CASM_clexmonte_run_covariance_functions

#include <limits>

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/state/Configuration.hh"
#include "casm/monte/results/ResultsAnalysisFunction.hh"
#include "casm/monte/state/StateSampler.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {

/// \brief Make variance analysis function (i.e. "heat_capacity")
monte::ResultsAnalysisFunction<Configuration> make_variance_f(
    std::string name, std::string description, std::string sampler_name,
    std::vector<std::string> component_names, std::vector<Index> shape,
    std::function<double(monte::Results<Configuration> const &)>
        make_normalization_constant_f);

/// \brief Make covariance analysis function (i.e. "mol_susc")
monte::ResultsAnalysisFunction<Configuration> make_covariance_f(
    std::string name, std::string description, std::string first_sampler_name,
    std::string second_sampler_name,
    std::vector<std::string> first_component_names,
    std::vector<std::string> second_component_names,
    std::function<double(monte::Results<Configuration> const &)>
        make_normalization_constant_f);

// /// \brief Make covariance of two components of sampled quantities
// ResultsAnalysisFunction<Configuration> make_component_covariance_f(...);

}  // namespace clexmonte
}  // namespace CASM

#endif
