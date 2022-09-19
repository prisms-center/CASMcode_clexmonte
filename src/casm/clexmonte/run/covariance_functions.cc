#include "casm/clexmonte/run/covariance_functions.hh"

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
///
/// \param name Name to give analysis function (ex. "heat_capacity")
/// \param description Description to give analysis function
/// \param sampler_name Name of state sampler function collecting
///     the data to calculate the variance of
/// \param component_names Names for components of the resulting
///     analysis
/// \param shape The shape of the resulting analysis. If the sampled
///     quantity is a scalar, then use `{}`, if the sampled quantity
///     is a vector, use `{v.size(), v.size()}`.
/// \param make_normalization_constant_f A function that returns the
///     the normalization constant. The analysis function will give
///    `result = var(sampler_name)/normalization_constant`.
///
monte::ResultsAnalysisFunction<Configuration> make_variance_f(
    std::string name, std::string description, std::string sampler_name,
    std::vector<std::string> component_names, std::vector<Index> shape,
    std::function<double(monte::Results<Configuration> const &)>
        make_normalization_constant_f) {
  return monte::ResultsAnalysisFunction<Configuration>(
      name, description, component_names, shape,
      [=](monte::Results<Configuration> const &results) {
        // validation of sampled data:
        auto it = results.samplers.find(sampler_name);
        if (it == results.samplers.end()) {
          std::stringstream msg;
          msg << "Results analysis error: " << name << " requires sampling "
              << sampler_name;
          throw std::runtime_error(msg.str());
        }
        auto const &sampler = *it->second;

        double normalization_constant = make_normalization_constant_f(results);

        monte::CompletionCheckResults const &completion_r =
            results.completion_check_results;
        monte::ConvergenceCheckResults const &convergence_r =
            completion_r.convergence_check_results;
        Index N_samples_for_statistics = convergence_r.N_samples_for_statistics;

        Index n = sampler.n_components();
        Eigen::VectorXd var = Eigen::VectorXd::Zero(n);
        for (Index i = 0; i < n; ++i) {
          if (N_samples_for_statistics == 0) {
            var(i) = std::numeric_limits<double>::quiet_NaN();
          } else {
            Eigen::VectorXd x =
                sampler.component(i).tail(N_samples_for_statistics);
            var(i) = monte::variance(x) / normalization_constant;
          }
        }

        return var;
      });
}

/// \brief Make covariance analysis function (i.e. "mol_susc")
///
/// \param name Name to give analysis function (ex. "mol_susc")
/// \param description Description to give analysis function
/// \param first_sampler_name Name of state sampler function
///     collecting the data to calculate the variance of.
/// \param second_sampler_name Name of state sampler function
///     collecting the data to calculate the variance of.
/// \param component_names Names for components of the resulting
///     analysis, with column-major unrolling
/// \param shape The shape of the resulting analysis. Use
///     `{first.size(), second.size()}`, with size=1 for scalar
///     quantities.
/// \param make_normalization_constant_f A function that returns the
///     the normalization constant. The analysis function will give
///    `result = cov(first, second)/normalization_constant`.
///
monte::ResultsAnalysisFunction<Configuration> make_covariance_f(
    std::string name, std::string description, std::string first_sampler_name,
    std::string second_sampler_name,
    std::vector<std::string> first_component_names,
    std::vector<std::string> second_component_names,
    std::function<double(monte::Results<Configuration> const &)>
        make_normalization_constant_f) {
  std::vector<std::string> cov_matrix_component_names;
  for (std::string col_name : second_component_names) {
    for (std::string row_name : first_component_names) {
      cov_matrix_component_names.push_back(row_name + "," + col_name);
    }
  }

  std::vector<Index> shape;
  shape.push_back(first_component_names.size());
  shape.push_back(second_component_names.size());

  return monte::ResultsAnalysisFunction<Configuration>(
      name, description, cov_matrix_component_names, shape,
      [=](monte::Results<Configuration> const &results) -> Eigen::VectorXd {
        // validation of sampled data:
        auto first_it = results.samplers.find(first_sampler_name);
        if (first_it == results.samplers.end()) {
          std::stringstream msg;
          msg << "Results analysis error: " << name << " requires sampling "
              << first_sampler_name;
          throw std::runtime_error(msg.str());
        }
        auto const &first_sampler = *first_it->second;

        // validation of sampled data:
        auto second_it = results.samplers.find(second_sampler_name);
        if (second_it == results.samplers.end()) {
          std::stringstream msg;
          msg << "Results analysis error: " << name << " requires sampling "
              << second_sampler_name;
          throw std::runtime_error(msg.str());
        }
        auto const &second_sampler = *second_it->second;

        double normalization_constant = make_normalization_constant_f(results);

        monte::CompletionCheckResults const &completion_r =
            results.completion_check_results;
        monte::ConvergenceCheckResults const &convergence_r =
            completion_r.convergence_check_results;
        Index N_samples_for_statistics = convergence_r.N_samples_for_statistics;

        Index m = first_sampler.n_components();
        Index n = second_sampler.n_components();
        Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(m, n);
        for (Index i = 0; i < m; ++i) {
          for (Index j = 0; j < n; ++j) {
            if (N_samples_for_statistics == 0) {
              cov(i, j) = std::numeric_limits<double>::quiet_NaN();
            } else {
              Eigen::VectorXd x =
                  first_sampler.component(i).tail(N_samples_for_statistics);
              Eigen::VectorXd y =
                  second_sampler.component(j).tail(N_samples_for_statistics);
              cov(i, j) = monte::covariance(x, y) / normalization_constant;
            }
          }
        }

        return monte::reshaped(cov);
      });
}

// /// \brief Make covariance of two components of sampled quantities
// template <typename ConfigType>
// ResultsAnalysisFunction<ConfigType> make_component_covariance_f(
//     SamplerComponent first, SamplerComponent second,
//     double normalization_constant, std::optional<std::string> name,
//     std::optional<std::string> description) {
//   if (!name.has_value()) {
//     std::stringstream ss;
//     ss << "cov(" << first.sampler_name << "(" << first.component_name << "),"
//        << second.sampler_name << "(" << second.component_name << "))";
//     name = ss.str();
//   }
//   if (!description.has_value()) {
//     std::stringstream ss;
//     ss << "Covariance of " << first.sampler_name << "(" <<
//     first.component_name
//        << ") and " << second.sampler_name << "(" << second.component_name
//        << ")";
//     description = ss.str();
//   }
//   return ResultsAnalysisFunction(
//       *name, *description, {},  // shape: empty vector for scalar
//       [=](Results<ConfigType> const &results) {
//         auto it_first = find_or_throw(results.samplers, first);
//         Eigen::VectorXd x =
//         it_first->second->component(first.component_index);
//
//         auto it_second = find_or_throw(results.samplers, second);
//         Eigen::VectorXd y =
//             it_second->second->component(second.component_index);
//
//         return covariance(x, y) / normalization_constant;
//       });
// }

}  // namespace clexmonte
}  // namespace CASM