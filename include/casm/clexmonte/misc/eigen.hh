#ifndef CASM_clexmonte_misc_eigen
#define CASM_clexmonte_misc_eigen

#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {
namespace clexmonte {

inline Eigen::VectorXd to_VectorXd(double value) {
  Eigen::VectorXd vec(1);
  vec(0) = value;
  return vec;
}

inline Eigen::VectorXd to_VectorXd(std::vector<double> const &value) {
  Eigen::VectorXd vec(value.size());
  for (int i = 0; i < value.size(); ++i) {
    vec(i) = value[i];
  }
  return vec;
}

}  // namespace clexmonte
}  // namespace CASM

#endif
