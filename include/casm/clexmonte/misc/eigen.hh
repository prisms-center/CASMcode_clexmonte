#ifndef CASM_clexmonte_misc_eigen
#define CASM_clexmonte_misc_eigen

#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {
namespace clexmonte {

inline Eigen::VectorXi to_VectorXi(int value) {
  return Eigen::VectorXi::Constant(1, value);
}

inline Eigen::VectorXi to_VectorXi(std::vector<int> const &value) {
  return Eigen::Map<const Eigen::VectorXi>(value.data(), value.size());
}

inline Eigen::VectorXl to_VectorXl(long value) {
  return Eigen::VectorXl::Constant(1, value);
}

inline Eigen::VectorXl to_VectorXl(std::vector<long> const &value) {
  return Eigen::Map<const Eigen::VectorXl>(value.data(), value.size());
}

inline Eigen::VectorXd to_VectorXd(double value) {
  return Eigen::VectorXd::Constant(1, value);
}

inline Eigen::VectorXd to_VectorXd(std::vector<double> const &value) {
  return Eigen::Map<const Eigen::VectorXd>(value.data(), value.size());
}

}  // namespace clexmonte
}  // namespace CASM

#endif
