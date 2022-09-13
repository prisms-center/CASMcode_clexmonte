#ifndef CASM_clexmonte_state_CorrMatchingPotential
#define CASM_clexmonte_state_CorrMatchingPotential

#include <optional>
#include <vector>

#include "casm/clexmonte/definitions.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}

namespace clexmonte {

struct CorrMatchingTarget {
  CorrMatchingTarget() : index(0), value(0.0), weight(1.0) {}

  CorrMatchingTarget(Index _index, double _value, double _weight)
      : index(_index), value(_value), weight(_weight) {}

  /// \brief Correlation index
  Index index;

  /// \brief Correlation target value
  double value;

  /// \brief Weight given to difference from target
  double weight;
};

/// \brief Parameters for a correlation-matching potential
///
/// Implements:
///
///     Epot = -w_{exact}*N_{exact} +
///         \sum_i v_i * | \Gamma_{j_i} - \Gamma^{target}_{j_i}) |,
///
/// where:
/// - N_{exact} is that maximum value such that
///   \Gamma_{j_i} - \Gamma^{target}_{j_i}) ~ 0 for all i < N_{exact}
/// - exact_matching_weight = w_{exact},
/// - targets[i].index = w_i
/// - targets[i].value = \Gamma^{target}_{j_i}
/// - targets[i].weight = v_i
struct CorrMatchingParams {
  CorrMatchingParams() : exact_matching_weight(0.0), tol(CASM::TOL) {}

  CorrMatchingParams(double _exact_matching_weight,
                     std::vector<CorrMatchingTarget> const &_targets,
                     double _tol)
      : exact_matching_weight(_exact_matching_weight),
        targets(_targets),
        tol(_tol) {}

  /// \brief Bias given for leading exactly matching correlations
  double exact_matching_weight;

  /// \brief Correlation matching targets
  std::vector<CorrMatchingTarget> targets;

  /// \brief Tolerance used to check for exactly matching correlations
  double tol;
};

double corr_matching_potential(Eigen::VectorXd const &corr,
                               CorrMatchingParams const &params);

Eigen::VectorXd make_corr_matching_error(Eigen::VectorXd const &corr,
                                         CorrMatchingParams const &params);

double delta_corr_matching_potential(Eigen::VectorXd const &corr,
                                     Eigen::VectorXd const &delta_corr,
                                     CorrMatchingParams const &params);

struct RandomAlloyCorrMatchingParams : public CorrMatchingParams {
  explicit RandomAlloyCorrMatchingParams(
      CorrCalculatorFunction _random_alloy_corr_f);

  RandomAlloyCorrMatchingParams(
      std::vector<Eigen::VectorXd> const &_sublattice_prob,
      CorrCalculatorFunction _random_alloy_corr_f,
      double _exact_matching_weight, std::vector<Index> _target_indices,
      std::optional<std::vector<double>> _target_weights = std::nullopt,
      double _tol = CASM::TOL);

  void update_targets();

  std::vector<Eigen::VectorXd> sublattice_prob;

  CorrCalculatorFunction random_alloy_corr_f;
};

// --- Used to convert to/from Eigen::VectorXd, for use by incremental
// conditions generator ---

Eigen::VectorXd to_VectorXd(CorrMatchingParams const &params);

Eigen::VectorXd to_VectorXd_increment(CorrMatchingParams const &params);

Eigen::VectorXd to_VectorXd(RandomAlloyCorrMatchingParams const &params);

Eigen::VectorXd to_VectorXd_increment(
    RandomAlloyCorrMatchingParams const &params);

/// \brief Specialize to help with conversions to/from VectorValueMap
template <typename ConstructedType>
struct ConditionsConstructor;

template <>
struct ConditionsConstructor<CorrMatchingParams> {
  static CorrMatchingParams from_VectorXd(Eigen::VectorXd const &v, double tol);
};

template <>
struct ConditionsConstructor<RandomAlloyCorrMatchingParams> {
  static RandomAlloyCorrMatchingParams from_VectorXd(
      Eigen::VectorXd const &v, xtal::BasicStructure const &prim,
      CorrCalculatorFunction random_alloy_corr_f, double tol);
};

// --- Print current values ---

void print_param(std::ostream &sout, std::string name,
                 std::optional<CorrMatchingParams> const &params);

void print_param(std::ostream &sout, std::string name,
                 std::optional<RandomAlloyCorrMatchingParams> const &params);

}  // namespace clexmonte
}  // namespace CASM

#endif
