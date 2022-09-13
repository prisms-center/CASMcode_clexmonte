#include "casm/clexmonte/state/CorrMatchingPotential.hh"

#include <iostream>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace clexmonte {

double corr_matching_potential(Eigen::VectorXd const &corr,
                               CorrMatchingParams const &params) {
  double Epot = 0;
  bool counting_n_exact = true;
  Index n_exact = 0;

  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential: target index out "
          "of range");
    }
    double value = corr(target.index);
    if (counting_n_exact) {
      if (CASM::almost_equal(value, target.value, params.tol)) {
        ++n_exact;
      } else {
        counting_n_exact = false;
      }
    }
    Epot += target.weight * std::abs(value - target.value);
  }
  Epot -= params.exact_matching_weight * n_exact;
  return Epot;
}

Eigen::VectorXd make_corr_matching_error(Eigen::VectorXd const &corr,
                                         CorrMatchingParams const &params) {
  Eigen::VectorXd corr_matching_error = Eigen::VectorXd::Zero(corr.size());

  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential: target index out "
          "of range");
    }
    corr_matching_error(target.index) = corr(target.index) - target.value;
  }
  return corr_matching_error;
}

double delta_corr_matching_potential(Eigen::VectorXd const &corr,
                                     Eigen::VectorXd const &delta_corr,
                                     CorrMatchingParams const &params) {
  if (corr.size() != delta_corr.size()) {
    throw std::runtime_error(
        "Error calculating correlations matching potential delta: corr and "
        "delta_corr size mismatch");
  }
  double dEpot = 0;
  bool counting_n_exact_1 = true;
  Index n_exact_1 = 0;
  bool counting_n_exact_2 = true;
  Index n_exact_2 = 0;
  for (auto const &target : params.targets) {
    if (target.index < 0 || target.index >= corr.size()) {
      throw std::runtime_error(
          "Error calculating correlations matching potential delta: target "
          "index out of range");
    }
    double value = corr(target.index);
    double dvalue = delta_corr(target.index);
    if (counting_n_exact_1) {
      if (CASM::almost_equal(value, target.value, params.tol)) {
        ++n_exact_1;
      } else {
        counting_n_exact_1 = false;
      }
    }
    if (counting_n_exact_2) {
      if (CASM::almost_equal(value + dvalue, target.value, params.tol)) {
        ++n_exact_2;
      } else {
        counting_n_exact_2 = false;
      }
    }
    dEpot += target.weight * (std::abs(value + dvalue - target.value) -
                              std::abs(value - target.value));
  }
  dEpot -= params.exact_matching_weight * (n_exact_2 - n_exact_1);
  return dEpot;
}

RandomAlloyCorrMatchingParams::RandomAlloyCorrMatchingParams(
    CorrCalculatorFunction _random_alloy_corr_f)
    : CorrMatchingParams(),
      sublattice_prob(),
      random_alloy_corr_f(_random_alloy_corr_f) {}

RandomAlloyCorrMatchingParams::RandomAlloyCorrMatchingParams(
    std::vector<Eigen::VectorXd> const &_sublattice_prob,
    CorrCalculatorFunction _random_alloy_corr_f, double _exact_matching_weight,
    std::vector<Index> _target_indices,
    std::optional<std::vector<double>> _target_weights, double _tol)
    : CorrMatchingParams(),
      sublattice_prob(_sublattice_prob),
      random_alloy_corr_f(_random_alloy_corr_f) {
  this->exact_matching_weight = _exact_matching_weight;
  this->tol = _tol;
  if (!_target_weights.has_value()) {
    _target_weights = std::vector<double>();
    for (Index i = 0; i < _target_indices.size(); ++i) {
      _target_weights->push_back(1.0);
    }
  }

  Index i = 0;
  for (auto index : _target_indices) {
    this->targets.emplace_back(index, 0.0, _target_weights->at(i));
    ++i;
  }
  this->update_targets();
}

void RandomAlloyCorrMatchingParams::update_targets() {
  Eigen::VectorXd random_alloy_corr =
      this->random_alloy_corr_f(sublattice_prob);
  for (auto &target : this->targets) {
    if (target.index >= random_alloy_corr.size()) {
      throw std::runtime_error(
          "Error in RandomAlloyCorrMatchingParams: correlation index out of "
          "range");
    }
    target.value = random_alloy_corr(target.index);
  }
}

/// \brief Convert CorrMatchingParams to an Eigen::VectorXd
///
/// Order:
/// - exact_matching_weight,
/// - targets[0].index, targets[0].value, targets[0].weight,
/// - targets[1].index, targets[1].value, targets[1].weight,
/// - ...
/// - targets[targets.size()-1].index, targets[targets.size()-1].value,
///   targets[targets.size()-1].weight
///
/// Note:
/// - `tol` is not saved to the vector
Eigen::VectorXd to_VectorXd(CorrMatchingParams const &params) {
  Eigen::VectorXd v(1 + 3 * params.targets.size());
  Index i = 0;
  v(i) = params.exact_matching_weight;
  ++i;
  for (auto const &target : params.targets) {
    v(i) = target.index;
    ++i;
    v(i) = target.value;
    ++i;
    v(i) = target.weight;
    ++i;
  }
  return v;
}

/// \brief Convert CorrMatchingParams to an Eigen::VectorXd that can be
///     added to increment CorrMatchingParams
///
/// Order (does not include targets[i].index):
/// - exact_matching_weight,
/// - 0, targets[0].value, targets[0].weight,
/// - 0, targets[1].value, targets[1].weight,
/// - ...
/// - 0, targets[targets.size()-1].value,
///   targets[targets.size()-1].weight
///
/// Note:
/// - `tol` is not saved to the vector
Eigen::VectorXd to_VectorXd_increment(CorrMatchingParams const &params) {
  Eigen::VectorXd v(1 + 3 * params.targets.size());
  Index i = 0;
  v(i) = params.exact_matching_weight;
  ++i;
  for (auto const &target : params.targets) {
    v(i) = 0.0;
    ++i;
    v(i) = target.value;
    ++i;
    v(i) = target.weight;
    ++i;
  }
  return v;
}

/// \brief Convert RandomAlloyCorrCorrMatchingParams to an Eigen::VectorXd
///
/// Order:
/// - exact_matching_weight,
/// - for each sublattice, b: sublattice_prob[b][0], sublattice_prob[b][1], ...,
/// sublattice_prob[b][sublattice_prob[b].size()-1]
/// - for each target, i: targets[i].index, targets[i].weight
///
/// Note:
/// - `tol` is not saved to the vector
Eigen::VectorXd to_VectorXd(RandomAlloyCorrMatchingParams const &params) {
  Index vsize = 1;
  for (auto const &sublat_prob : params.sublattice_prob) {
    vsize += sublat_prob.size();
  }
  vsize += params.targets.size();

  Eigen::VectorXd v(vsize);
  Index i = 0;
  v(i) = params.exact_matching_weight;
  ++i;
  for (auto const &sublat_prob : params.sublattice_prob) {
    for (auto const &value : sublat_prob) {
      v(i) = value;
      ++i;
    }
  }
  for (auto const &target : params.targets) {
    v(i) = target.index;
    ++i;
    v(i) = target.weight;
    ++i;
  }
  return v;
}

/// \brief Convert RandomAlloyCorrCorrMatchingParams to an Eigen::VectorXd
///
/// Order:
/// - exact_matching_weight,
/// - for each sublattice, b: sublattice_prob[b][0], sublattice_prob[b][1], ...,
/// sublattice_prob[b][sublattice_prob[b].size()-1]
/// - for each target, i: 0, 0
///
/// Note:
/// - `tol` is not saved to the vector
Eigen::VectorXd to_VectorXd_increment(
    RandomAlloyCorrMatchingParams const &params) {
  Index vsize = 1;
  for (auto const &sublat_prob : params.sublattice_prob) {
    vsize += sublat_prob.size();
  }
  vsize += params.targets.size();

  Eigen::VectorXd v(vsize);
  Index i = 0;
  v(i) = params.exact_matching_weight;
  ++i;
  for (auto const &sublat_prob : params.sublattice_prob) {
    for (auto const &value : sublat_prob) {
      v(i) = value;
      ++i;
    }
  }
  for (Index t = 0; t < params.targets.size(); ++t) {
    v(i) = 0;
    ++i;
    v(i) = params.targets[t].weight;
    ++i;
  }
  return v;
}

/// \brief Convert an Eigen::VectorXd to CorrMatchingParams
///
/// Inverse process as `to_VectorXd(CorrMatchingParams const &params)`
///
/// Note:
/// - `tol` must be provided separately
CorrMatchingParams ConditionsConstructor<CorrMatchingParams>::from_VectorXd(
    Eigen::VectorXd const &v, double tol) {
  double exact_matching_weight = 0;
  if (v.size() < 1) {
    throw std::runtime_error(
        "Error reading CorrMatchingParams: incorrect vector size, size = 0");
  }
  Index i = 0;
  exact_matching_weight = v(0);
  ++i;

  std::vector<CorrMatchingTarget> targets;
  while (i != v.size()) {
    if (i + 3 > v.size()) {
      throw std::runtime_error(
          "Error reading CorrMatchingParams: incorrect vector size, size != 1  "
          "+ 3*n_targets");
    }
    Index index = std::lround(v(i));
    ++i;
    double value = v(i);
    ++i;
    double weight = v(i);
    ++i;
    targets.emplace_back(index, value, weight);
  }

  return CorrMatchingParams{exact_matching_weight, targets, tol};
}

/// \brief Convert an Eigen::VectorXd to RandomAlloyCorrMatchingParams
///
/// Inverse process as `to_VectorXd(RandomAlloyCorrMatchingParams const
/// &params)`
///
/// Note:
/// - `tol` must be provided separately
RandomAlloyCorrMatchingParams
ConditionsConstructor<RandomAlloyCorrMatchingParams>::from_VectorXd(
    Eigen::VectorXd const &v, xtal::BasicStructure const &prim,
    CorrCalculatorFunction random_alloy_corr_f, double tol) {
  double exact_matching_weight = 0;
  if (v.size() < 1) {
    throw std::runtime_error(
        "Error reading RandomAlloyCorrMatchingParams: incorrect vector size, "
        "size = 0");
  }
  Index i = 0;
  exact_matching_weight = v(0);
  ++i;

  std::vector<Eigen::VectorXd> sublattice_prob;
  for (auto const &site : prim.basis()) {
    Eigen::VectorXd v_sublat(site.occupant_dof().size());
    for (Index j = 0; j < site.occupant_dof().size(); ++j) {
      if (i == v.size()) {
        throw std::runtime_error(
            "Error reading RandomAlloyCorrMatchingParams: incorrect vector "
            "size");
      }
      v_sublat(j) = v(i);
      ++i;
    }
    sublattice_prob.push_back(v_sublat);
  }

  std::vector<Index> target_indices;
  std::vector<double> target_weights;
  while (i < v.size()) {
    target_indices.push_back(v(i));
    ++i;
    target_weights.push_back(v(i));
    ++i;
  }

  return RandomAlloyCorrMatchingParams{
      sublattice_prob, random_alloy_corr_f, exact_matching_weight,
      target_indices,  target_weights,      tol};
}

void print_param(std::ostream &sout, std::string name,
                 std::optional<CorrMatchingParams> const &params) {
  if (!params.has_value()) {
    return;
  }
  sout << name << ".exact_matching_weight: " << params->exact_matching_weight
       << "\n";
  for (auto const &target : params->targets) {
    sout << name << ".target: {index=" << target.index
         << ", value=" << target.value << ", weight=" << target.weight << "}\n";
  }
}

void print_param(std::ostream &sout, std::string name,
                 std::optional<RandomAlloyCorrMatchingParams> const &params) {
  if (!params.has_value()) {
    return;
  }
  sout << name << ".sublattice_prob: " << std::endl;
  Index b = 0;
  for (auto const &prob : params->sublattice_prob) {
    sout << "  " << b << ": " << prob.transpose() << std::endl;
    ++b;
  }
  sout << name << ".exact_matching_weight: " << params->exact_matching_weight
       << "\n";
}

}  // namespace clexmonte
}  // namespace CASM
