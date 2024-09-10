#ifndef CASM_clexmonte_diffusion_calculations
#define CASM_clexmonte_diffusion_calculations

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace clexmonte {

/// \brief Handle unrolling convention for L_{ij}
struct CollectiveIsotropicCounter {
  std::vector<std::string> name_list;
  Index i;
  Index j;

  CollectiveIsotropicCounter(std::vector<std::string> _name_list)
      : name_list(_name_list) {
    reset();
  }

  void reset() {
    i = 0;
    j = 0;
  }

  std::string component_name() const {
    return name_list[i] + "," + name_list[j];
  }

  void advance() {
    ++j;
    if (j == name_list.size()) {
      ++i;
      j = i;
    }
  }

  bool is_valid() const { return i != name_list.size(); }
};

/// \brief Handle unrolling convention for D^{*}_{i}
struct IndividualIsotropicCounter {
  std::vector<std::string> name_list;
  Index i;

  IndividualIsotropicCounter(std::vector<std::string> _name_list)
      : name_list(_name_list) {
    reset();
  }

  void reset() { i = 0; }

  std::string component_name() const { return name_list[i]; }

  void advance() { ++i; }

  bool is_valid() const { return i != name_list.size(); }
};

/// \brief Base class to handle unrolling convention for directions
struct AnisotropicCounter {
  std::vector<std::string> name_list;
  std::vector<std::string> dirs;
  std::vector<Index> index_alpha;
  std::vector<Index> index_beta;

  Index dir_index;
  Index alpha;
  Index beta;

  AnisotropicCounter(std::vector<std::string> _name_list)
      : name_list(_name_list),
        dirs({"x", "y", "z"}),
        index_alpha({0, 1, 2, 1, 0, 0}),
        index_beta({0, 1, 2, 2, 2, 1}) {
    dir_reset();
  }

  void dir_reset() {
    dir_index = 0;
    alpha = index_alpha[dir_index];
    beta = index_beta[dir_index];
  }

  std::string dir_component_name() const {
    return dirs[alpha] + "," + dirs[beta];
  }
};

/// \brief Handle unrolling convention for L_{ij\alpha\beta}
struct CollectiveAnisotropicCounter : public AnisotropicCounter {
  Index i;
  Index j;

  CollectiveAnisotropicCounter(std::vector<std::string> _name_list)
      : AnisotropicCounter(_name_list) {
    reset();
  }

  void reset() {
    i = 0;
    j = 0;
    this->dir_reset();
  }

  std::string component_name() const {
    return name_list[i] + "," + name_list[j] + "," + dir_component_name();
  }

  void advance() {
    ++dir_index;
    if (dir_index == 6) {
      ++j;
      dir_index = 0;
    }
    if (j == name_list.size()) {
      ++i;
      j = i;
      dir_index = 0;
    }
    alpha = index_alpha[dir_index];
    beta = index_beta[dir_index];
  }

  bool is_valid() const { return i != name_list.size(); }
};

/// \brief Handle unrolling convention for D^{*}_{i\alpha\beta}
struct IndividualAnisotropicCounter : public AnisotropicCounter {
  Index i;

  IndividualAnisotropicCounter(std::vector<std::string> _name_list)
      : AnisotropicCounter(_name_list) {
    reset();
  }

  void reset() {
    i = 0;
    this->dir_reset();
  }

  std::string component_name() const {
    return name_list[i] + "," + dir_component_name();
  }

  void advance() {
    ++dir_index;
    if (dir_index == 6) {
      ++i;
      dir_index = 0;
    }
    alpha = index_alpha[dir_index];
    beta = index_beta[dir_index];
  }

  bool is_valid() const { return i != name_list.size(); }
};

/// \brief Get unrolled component names
template <typename CounterType>
static std::vector<std::string> make_component_names(
    std::vector<std::string> name_list) {
  CounterType counter(name_list);
  std::vector<std::string> component_names;
  while (counter.is_valid()) {
    component_names.push_back(counter.component_name());
    counter.advance();
  }
  return component_names;
}

/// \brief Returns the mean collective atom isotropic squared displacement
///
/// This is:
///
///     \frac{1}{N} (\sum_\zeta \Delta R^\zeta_{i}) \dot (\sum_\zeta \Delta
///     R^\zeta_{j})
///
/// where i is atom type, \zeta is index of atoms of type i, and N is the total
/// number of atoms.
///
/// \returns A vector, v, containing values unrolled in the order specified by
///      `make_components_names<CollectiveIsotropicCounter>(atom_name_list)`.
///
/// Notes:
/// - The unrolled values order should be:
///     (00, 01, 02, ..., 10, ..., nn, ...)
///
inline Eigen::VectorXd mean_R_squared_collective_isotropic(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R) {
  std::vector<Eigen::Vector3d> sumR(atom_name_list.size(),
                                    Eigen::Vector3d::Zero());
  for (Index atom_index = 0; atom_index < delta_R.cols(); ++atom_index) {
    sumR[atom_name_index_list[atom_index]] += delta_R.col(atom_index);
  }

  CollectiveIsotropicCounter counter(atom_name_list);
  std::vector<double> v;
  double N = delta_R.cols();
  while (counter.is_valid()) {
    v.push_back(sumR[counter.i].dot(sumR[counter.j]) / N);
    counter.advance();
  }
  return to_VectorXd(v);
}

inline Eigen::VectorXd L_isotropic_sample(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R,
    double dim, double delta_time, double volume, double temperature) {
  // double n_unitcells =
  // get_transformation_matrix_to_super(state).determinant(); double volume
  // = n_unitcells * get_prim_basicstructure(system)->lattice().volume();
  double n_atoms = atom_name_index_list.size();
  Eigen::VectorXd RiRj = mean_R_squared_collective_isotropic(
                             atom_name_list, atom_name_index_list, delta_R) *
                         n_atoms;
  return RiRj / (2 * dim * delta_time * volume * CASM::KB * temperature);
}

/// \brief Returns the mean collective atom anisotropic squared displacement
///
/// This is:
///
///     \frac{1}{N} (\sum_\zeta \Delta R^\zeta_{i,\alpha}) * (\sum_\zeta \Delta
///     R^\zeta_{j,\beta})
///
/// where i is atom type, \zeta is index of atoms of type i, and \alpha, \beta
/// indicate directions (x, y, z), and N is the total number of atoms.
///
/// \returns A vector, v, containing values unrolled in the order specified by
///      `make_components_names<CollectiveAnisotropicCounter>(atom_name_list)`.
///
/// Notes:
/// - The unrolled values order should be:
///     (00xx, 00yy, 00zz, 00yz, 00xz, 00yx, 01xx, ..., 02xx, ..., 10xx, ...,
///     nnxx, ...)
///
inline Eigen::VectorXd mean_R_squared_collective_anisotropic(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R) {
  std::vector<Eigen::Vector3d> sumR(atom_name_list.size(),
                                    Eigen::Vector3d::Zero());
  for (Index atom_index = 0; atom_index < delta_R.cols(); ++atom_index) {
    sumR[atom_name_index_list[atom_index]] += delta_R.col(atom_index);
  }

  CollectiveAnisotropicCounter counter(atom_name_list);
  std::vector<double> v;
  double N = delta_R.cols();
  while (counter.is_valid()) {
    v.push_back(sumR[counter.i](counter.alpha) * sumR[counter.j](counter.beta) /
                N);
    counter.advance();
  }
  return to_VectorXd(v);
}

inline Eigen::VectorXd L_anisotropic_sample(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R,
    double delta_time, double volume, double temperature) {
  // double n_unitcells =
  // get_transformation_matrix_to_super(state).determinant(); double volume
  // = n_unitcells * get_prim_basicstructure(system)->lattice().volume();
  double n_atoms = atom_name_index_list.size();
  Eigen::VectorXd RiaRjb = mean_R_squared_collective_anisotropic(
                               atom_name_list, atom_name_index_list, delta_R) *
                           n_atoms;
  return RiaRjb / (2 * delta_time * volume * CASM::KB * temperature);
}

/// \brief Returns the mean individual atom isotropic squared displacement
///
/// This is:
///
///     \frac{1}{N_i} \sum_\zeta \Delta R^\zeta_{i} \dot \Delta R^\zeta_{i}
///
/// where i is atom type, \zeta is an index for atoms of type i, and N_i is
/// the number of atoms of type i.
///
/// \returns A vector, v, containing values unrolled in the order specified by
///      `make_components_names<IndividualIsotropicCounter>(atom_name_list)`.
///
/// Notes:
/// - The unrolled values order should be:
///     (0xx, 0yy, 0zz, 0yz, 0xz, 0yx, 1xx, ..., 2xx, ..., ...)
///
inline Eigen::VectorXd mean_R_squared_individual_isotropic(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R) {
  std::vector<double> sumRR(atom_name_list.size(), 0.0);
  Eigen::VectorXd N = Eigen::VectorXd::Zero(atom_name_list.size());
  for (Index atom_index = 0; atom_index < delta_R.cols(); ++atom_index) {
    sumRR[atom_name_index_list[atom_index]] +=
        delta_R.col(atom_index).dot(delta_R.col(atom_index));
    N(atom_name_index_list[atom_index]) += 1.0;
  }

  IndividualIsotropicCounter counter(atom_name_list);
  std::vector<double> v;
  while (counter.is_valid()) {
    v.push_back(sumRR[counter.i] / N(counter.i));
    counter.advance();
  }
  return to_VectorXd(v);
}

inline Eigen::VectorXd D_tracer_isotropic_sample(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R,
    double dim, double delta_time) {
  Eigen::VectorXd RiRi = mean_R_squared_individual_isotropic(
      atom_name_list, atom_name_index_list, delta_R);
  return RiRi / (2 * dim * delta_time);
}

/// \brief Returns the mean individual atom anisotropic squared displacement
///
/// This is:
///
///     \frac{1}{N_i} \sum_\zeta \Delta R^\zeta_{i,\alpha} * \Delta
///     R^\zeta_{i,\beta}
///
/// where i is atom type, \zeta is an index for atoms of type i, and \alpha,
/// \beta indicate directions (x, y, z), and N_i is the number of atoms of type
/// i.
///
/// \returns A vector, v, containing values unrolled in the order specified by
///      `make_components_names<IndividualAnisotropicCounter>(atom_name_list)`.
///
/// Notes:
/// - The unrolled values order should be:
///     (0xx, 0yy, 0zz, 0yz, 0xz, 0yx, 1xx, ..., 2xx, ..., ...)
///
inline Eigen::VectorXd mean_R_squared_individual_anisotropic(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R) {
  std::vector<Eigen::Matrix3d> sumRR(atom_name_list.size(),
                                     Eigen::Matrix3d::Zero());
  Eigen::VectorXd N = Eigen::VectorXd::Zero(atom_name_list.size());
  for (Index atom_index = 0; atom_index < delta_R.cols(); ++atom_index) {
    sumRR[atom_name_index_list[atom_index]] +=
        delta_R.col(atom_index) * delta_R.col(atom_index).transpose();
    N(atom_name_index_list[atom_index]) += 1.0;
  }

  IndividualAnisotropicCounter counter(atom_name_list);
  std::vector<double> v;
  while (counter.is_valid()) {
    v.push_back(sumRR[counter.i](counter.alpha, counter.beta) / N(counter.i));
    counter.advance();
  }
  return to_VectorXd(v);
}

inline Eigen::VectorXd D_tracer_anisotropic_sample(
    std::vector<std::string> atom_name_list,
    std::vector<Index> atom_name_index_list, Eigen::MatrixXd const &delta_R,
    double delta_time) {
  Eigen::VectorXd RiaRib = mean_R_squared_individual_anisotropic(
      atom_name_list, atom_name_index_list, delta_R);
  return RiaRib / (2 * delta_time);
}

}  // namespace clexmonte
}  // namespace CASM

#endif
