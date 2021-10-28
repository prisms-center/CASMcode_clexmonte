#ifndef CASM_clexmonte_misc_Matrix3lCompare
#define CASM_clexmonte_misc_Matrix3lCompare

#include <algorithm>

#include "casm/global/eigen.hh"

namespace CASM {
namespace clexmonte {

struct Matrix3lCompare {
  bool operator()(Eigen::Matrix3l const &lhs,
                  Eigen::Matrix3l const &rhs) const {
    return std::lexicographical_compare(lhs.data(), lhs.data() + lhs.size(),
                                        rhs.data(), rhs.data() + rhs.size());
  }
};

}  // namespace clexmonte
}  // namespace CASM

#endif
