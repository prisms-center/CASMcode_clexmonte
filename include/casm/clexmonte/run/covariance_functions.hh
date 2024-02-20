#ifndef CASM_clexmonte_run_covariance_functions
#define CASM_clexmonte_run_covariance_functions

#include "casm/clexmonte/definitions.hh"

namespace CASM {
namespace clexmonte {

/// \brief Make variance analysis function (i.e. "heat_capacity")
results_analysis_function_type make_variance_f(
    std::string name, std::string description, std::string sampler_name,
    std::vector<std::string> component_names, std::vector<Index> shape,
    std::function<double()> make_normalization_constant_f);

/// \brief Make covariance analysis function (i.e. "mol_susc")
results_analysis_function_type make_covariance_f(
    std::string name, std::string description, std::string first_sampler_name,
    std::string second_sampler_name,
    std::vector<std::string> first_component_names,
    std::vector<std::string> second_component_names,
    std::function<double()> make_normalization_constant_f);

}  // namespace clexmonte
}  // namespace CASM

#endif
