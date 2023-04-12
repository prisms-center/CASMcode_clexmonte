#include <cstddef>

#include "casm/clexulator/BaseClexulator.hh"
#include "casm/clexulator/BasicClexParamPack.hh"
#include "casm/global/eigen.hh"

/****** PROJECT SPECIFICATIONS ******

         ****** prim.json ******

{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],
      "occupants" : [ "A", "B", "Va" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "lattice_vectors" : [
    [ 0.000000000000, 2.000000000000, 2.000000000000 ],
    [ 2.000000000000, 0.000000000000, 2.000000000000 ],
    [ 2.000000000000, 2.000000000000, 0.000000000000 ]
  ],
  "title" : "FCC_binary_vacancy"
}

        ****** bspecs.json ******

{
  "basis_function_specs" : {
    "dof_specs" : {
      "occ" : {
        "site_basis_functions" : "OCCUPATION"
      }
    },
    "dofs" : [ "occ" ],
    "global_max_poly_order" : -1,
    "param_pack_type" : "DEFAULT"
  },
  "cluster_specs" : {
    "method" : "periodic_max_length",
    "params" : {
      "generating_group" : [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47 ], "orbit_branch_specs" : {
        "0" : {
          "max_length" : 0.000000000000
        },
        "1" : {
          "max_length" : 0.000000000000
        },
        "2" : {
          "max_length" : 4.100000000000
        }
      }
    }
  }
}

**/

/// \brief Returns a clexulator::BaseClexulator* owning a
/// FCC_binary_vacancy_Clexulator_default
extern "C" CASM::clexulator::BaseClexulator *
make_FCC_binary_vacancy_Clexulator_default();

namespace CASM {
namespace clexulator {

/****** GENERATED CLEXPARAMPACK DEFINITION ******/

typedef BasicClexParamPack ParamPack;

/****** GENERATED CLEXULATOR DEFINITION ******/

class FCC_binary_vacancy_Clexulator_default
    : public clexulator::BaseClexulator {
 public:
  FCC_binary_vacancy_Clexulator_default();

  ~FCC_binary_vacancy_Clexulator_default();

  ClexParamPack const &param_pack() const override { return m_params; }

  ClexParamPack &param_pack() override { return m_params; }

  template <typename Scalar>
  Scalar eval_bfunc_0_0() const;

  template <typename Scalar>
  Scalar eval_bfunc_1_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_1_1() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_1_0_at_0() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_1_1_at_0() const;

  template <typename Scalar>
  Scalar site_deval_bfunc_1_0_at_0(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_1_1_at_0(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar eval_bfunc_2_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_2_1() const;
  template <typename Scalar>
  Scalar eval_bfunc_2_2() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_2_0_at_0() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_1_at_0() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_2_at_0() const;

  template <typename Scalar>
  Scalar site_deval_bfunc_2_0_at_0(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_1_at_0(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_2_at_0(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar eval_bfunc_3_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_3_1() const;
  template <typename Scalar>
  Scalar eval_bfunc_3_2() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_3_0_at_0() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_1_at_0() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_2_at_0() const;

  template <typename Scalar>
  Scalar site_deval_bfunc_3_0_at_0(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_1_at_0(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_2_at_0(int occ_i, int occ_f) const;

 private:
  // ParamPack object, which stores temporary data for calculations
  mutable ParamPack m_params;

  // typedef for method pointers of scalar type double
  typedef double (
      FCC_binary_vacancy_Clexulator_default::*BasisFuncPtr_0)() const;

  // typedef for method pointers
  typedef double (FCC_binary_vacancy_Clexulator_default::*DeltaBasisFuncPtr_0)(
      int, int) const;

  // array of pointers to member functions for calculating basis functions of
  // scalar type double
  BasisFuncPtr_0 m_orbit_func_table_0[9];

  // array of pointers to member functions for calculating flower functions of
  // scalar type double
  BasisFuncPtr_0 m_flower_func_table_0[1][9];

  // array of pointers to member functions for calculating DELTA flower
  // functions of scalar type double
  DeltaBasisFuncPtr_0 m_delta_func_table_0[1][9];

  // Occupation Function tables for basis sites in asymmetric unit 0:
  //   - basis site 0:
  double m_occ_func_0_0[3];
  double m_occ_func_0_1[3];

  // ClexParamPack allocation for evaluated correlations
  ParamPack::Key m_corr_param_key;
  // ClexParamPack allocation for DoF occ
  ParamPack::Key m_occ_site_func_param_key;

  /// \brief Clone the FCC_binary_vacancy_Clexulator_default
  BaseClexulator *_clone() const override {
    return new FCC_binary_vacancy_Clexulator_default(*this);
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  /// Result is recorded in ClexParamPack
  void _calc_global_corr_contribution() const override;

  /// \brief Calculate contribution to global correlations from one unit cell
  /// /// Result is recorded in double array starting at corr_begin
  void _calc_global_corr_contribution(double *corr_begin) const override;

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell into ClexParamPack Result is recorded in ClexParamPack
  void _calc_restricted_global_corr_contribution(
      size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate contribution to select global correlations from one unit
  /// cell Result is recorded in double array starting at corr_begin
  void _calc_restricted_global_corr_contribution(
      double *corr_begin, size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate point correlations about neighbor site 'nlist_ind'
  /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
  /// For local clexulators, 'nlist_ind' ranges over all sites in the
  /// neighborhood Result is recorded in ClexParamPack
  void _calc_point_corr(int nlist_ind) const override;

  /// \brief Calculate point correlations about neighbor site 'nlist_ind'
  /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
  /// For local clexulators, 'nlist_ind' ranges over all sites in the
  /// neighborhood Result is recorded in double array starting at corr_begin
  void _calc_point_corr(int nlist_ind, double *corr_begin) const override;

  /// \brief Calculate select point correlations about neighbor site 'nlist_ind'
  /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
  /// For local clexulators, 'nlist_ind' ranges over all sites in the
  /// neighborhood Result is recorded in ClexParamPack
  void _calc_restricted_point_corr(
      int nlist_ind, size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate select point correlations about neighbor site 'nlist_ind'
  /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
  /// For local clexulators, 'nlist_ind' ranges over all sites in the
  /// neighborhood Result is recorded in double array starting at corr_begin
  void _calc_restricted_point_corr(
      int nlist_ind, double *corr_begin, size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant at neighbor site 'nlist_ind' For global clexulators, 'nlist_ind'
  /// only ranges over sites in the cell For local clexulators, 'nlist_ind'
  /// ranges over all sites in the neighborhood Result is recorded in
  /// ClexParamPack
  void _calc_delta_point_corr(int nlist_ind, int occ_i,
                              int occ_f) const override;

  /// \brief Calculate the change in point correlations due to changing an
  /// occupant at neighbor site 'nlist_ind' For global clexulators, 'nlist_ind'
  /// only ranges over sites in the cell For local clexulators, 'nlist_ind'
  /// ranges over all sites in the neighborhood Result is recorded in double
  /// array starting at corr_begin
  void _calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f,
                              double *corr_begin) const override;

  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant at neighbor site 'nlist_ind' For global clexulators,
  /// 'nlist_ind' only ranges over sites in the cell For local clexulators,
  /// 'nlist_ind' ranges over all sites in the neighborhood Result is recorded
  /// in ClexParamPack
  void _calc_restricted_delta_point_corr(
      int nlist_ind, int occ_i, int occ_f, size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  /// \brief Calculate the change in select point correlations due to changing
  /// an occupant at neighbor site 'nlist_ind' For global clexulators,
  /// 'nlist_ind' only ranges over sites in the cell For local clexulators,
  /// 'nlist_ind' ranges over all sites in the neighborhood Result is recorded
  /// in double array starting at corr_begin
  void _calc_restricted_delta_point_corr(
      int nlist_ind, int occ_i, int occ_f, double *corr_begin,
      size_type const *ind_list_begin,
      size_type const *ind_list_end) const override;

  template <typename Scalar>
  void _global_prepare() const;

  template <typename Scalar>
  void _point_prepare(int nlist_ind) const;

  // Occupation Function evaluators and accessors for basis site 0:
  double const &eval_occ_func_0_0(const int &nlist_ind) const {
    return m_occ_func_0_0[_occ(nlist_ind)];
  }

  double const &occ_func_0_0(const int &nlist_ind) const {
    return m_params.read(m_occ_site_func_param_key, 0, nlist_ind);
  }
  double const &eval_occ_func_0_1(const int &nlist_ind) const {
    return m_occ_func_0_1[_occ(nlist_ind)];
  }

  double const &occ_func_0_1(const int &nlist_ind) const {
    return m_params.read(m_occ_site_func_param_key, 1, nlist_ind);
  }

  // default functions for basis function evaluation
  template <typename Scalar>
  Scalar zero_func() const {
    return Scalar(0.0);
  }

  template <typename Scalar>
  Scalar zero_func(int, int) const {
    return Scalar(0.0);
  }
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FCC_binary_vacancy_Clexulator_default::FCC_binary_vacancy_Clexulator_default()
    : BaseClexulator(19, 9, 1) {
  m_occ_func_0_0[0] = 0.0000000000, m_occ_func_0_0[1] = 1.0000000000,
  m_occ_func_0_0[2] = 0.0000000000;

  m_occ_func_0_1[0] = 0.0000000000, m_occ_func_0_1[1] = 0.0000000000,
  m_occ_func_0_1[2] = 1.0000000000;

  m_occ_site_func_param_key = m_params.allocate("occ_site_func", 2, 19, true);

  m_corr_param_key = m_params.allocate("corr", corr_size(), 1, false);

  m_orbit_func_table_0[0] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_0_0<double>;
  m_orbit_func_table_0[1] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_1_0<double>;
  m_orbit_func_table_0[2] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_1_1<double>;
  m_orbit_func_table_0[3] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_0<double>;
  m_orbit_func_table_0[4] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_1<double>;
  m_orbit_func_table_0[5] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_2<double>;
  m_orbit_func_table_0[6] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_0<double>;
  m_orbit_func_table_0[7] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_1<double>;
  m_orbit_func_table_0[8] =
      &FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_2<double>;

  m_flower_func_table_0[0][0] =
      &FCC_binary_vacancy_Clexulator_default::zero_func<double>;
  m_flower_func_table_0[0][1] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_1_0_at_0<double>;
  m_flower_func_table_0[0][2] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_1_1_at_0<double>;
  m_flower_func_table_0[0][3] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_0_at_0<double>;
  m_flower_func_table_0[0][4] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_1_at_0<double>;
  m_flower_func_table_0[0][5] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_2_at_0<double>;
  m_flower_func_table_0[0][6] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_0_at_0<double>;
  m_flower_func_table_0[0][7] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_1_at_0<double>;
  m_flower_func_table_0[0][8] =
      &FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_2_at_0<double>;

  m_delta_func_table_0[0][0] =
      &FCC_binary_vacancy_Clexulator_default::zero_func<double>;
  m_delta_func_table_0[0][1] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_1_0_at_0<double>;
  m_delta_func_table_0[0][2] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_1_1_at_0<double>;
  m_delta_func_table_0[0][3] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_0_at_0<double>;
  m_delta_func_table_0[0][4] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_1_at_0<double>;
  m_delta_func_table_0[0][5] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_2_at_0<double>;
  m_delta_func_table_0[0][6] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_0_at_0<double>;
  m_delta_func_table_0[0][7] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_1_at_0<double>;
  m_delta_func_table_0[0][8] =
      &FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_2_at_0<double>;

  m_weight_matrix.row(0) << 2, 1, 1;
  m_weight_matrix.row(1) << 1, 2, 1;
  m_weight_matrix.row(2) << 1, 1, 2;

  m_sublat_indices = std::set<int>{0};

  m_n_sublattices = 1;

  m_neighborhood = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, -1, 1), xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),  xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),  xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, -1, 0),  xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),  xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),   xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),   xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),  xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -1),  xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, -1)};

  m_orbit_neighborhood.resize(corr_size());
  m_orbit_site_neighborhood.resize(corr_size());
  m_orbit_neighborhood[1] = std::set<xtal::UnitCell>{xtal::UnitCell(0, 0, 0)};
  m_orbit_neighborhood[2] = m_orbit_neighborhood[1];

  m_orbit_site_neighborhood[1] =
      std::set<xtal::UnitCellCoord>{xtal::UnitCellCoord(0, 0, 0, 0)};
  m_orbit_site_neighborhood[2] = m_orbit_site_neighborhood[1];

  m_orbit_neighborhood[3] = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, 0, 0), xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, 0), xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1), xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),  xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1), xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0), xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0)};
  m_orbit_neighborhood[4] = m_orbit_neighborhood[3];
  m_orbit_neighborhood[5] = m_orbit_neighborhood[3];

  m_orbit_site_neighborhood[3] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, -1, 0, 0), xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, 0), xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1), xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),  xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1), xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0), xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0)};
  m_orbit_site_neighborhood[4] = m_orbit_site_neighborhood[3];
  m_orbit_site_neighborhood[5] = m_orbit_site_neighborhood[3];

  m_orbit_neighborhood[6] = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, -1, 1), xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 1),  xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -1, -1), xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 1, -1)};
  m_orbit_neighborhood[7] = m_orbit_neighborhood[6];
  m_orbit_neighborhood[8] = m_orbit_neighborhood[6];

  m_orbit_site_neighborhood[6] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, -1, -1, 1), xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 1),  xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -1, -1), xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 1, -1)};
  m_orbit_site_neighborhood[7] = m_orbit_site_neighborhood[6];
  m_orbit_site_neighborhood[8] = m_orbit_site_neighborhood[6];
}

FCC_binary_vacancy_Clexulator_default::
    ~FCC_binary_vacancy_Clexulator_default() {
  // nothing here for now
}

/// \brief Calculate contribution to global correlations from one unit cell
void FCC_binary_vacancy_Clexulator_default::_calc_global_corr_contribution(
    double *corr_begin) const {
  _calc_global_corr_contribution();
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate contribution to global correlations from one unit cell
void FCC_binary_vacancy_Clexulator_default::_calc_global_corr_contribution()
    const {
  m_params.pre_eval();
  {
    _global_prepare<double>();
    for (size_type i = 0; i < corr_size(); i++) {
      ParamPack::Val<double>::set(m_params, m_corr_param_key, i,
                                  (this->*m_orbit_func_table_0[i])());
    }
  }
  m_params.post_eval();
}

/// \brief Calculate contribution to select global correlations from one unit
/// cell
void FCC_binary_vacancy_Clexulator_default::
    _calc_restricted_global_corr_contribution(
        double *corr_begin, size_type const *ind_list_begin,
        size_type const *ind_list_end) const {
  _calc_restricted_global_corr_contribution(ind_list_begin, ind_list_end);
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(
        m_params, m_corr_param_key, *ind_list_begin);
  }
}

/// \brief Calculate contribution to select global correlations from one unit
/// cell
void FCC_binary_vacancy_Clexulator_default::
    _calc_restricted_global_corr_contribution(
        size_type const *ind_list_begin, size_type const *ind_list_end) const {
  m_params.pre_eval();
  {
    _global_prepare<double>();
    for (; ind_list_begin < ind_list_end; ind_list_begin++) {
      ParamPack::Val<double>::set(
          m_params, m_corr_param_key, *ind_list_begin,
          (this->*m_orbit_func_table_0[*ind_list_begin])());
    }
  }
  m_params.post_eval();
}

/// \brief Calculate point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_default::_calc_point_corr(
    int nlist_ind, double *corr_begin) const {
  _calc_point_corr(nlist_ind);
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_default::_calc_point_corr(
    int nlist_ind) const {
  m_params.pre_eval();
  {
    _point_prepare<double>(nlist_ind);
    for (size_type i = 0; i < corr_size(); i++) {
      ParamPack::Val<double>::set(
          m_params, m_corr_param_key, i,
          (this->*m_flower_func_table_0[nlist_ind][i])());
    }
  }
  m_params.post_eval();
}

/// \brief Calculate select point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_default::_calc_restricted_point_corr(
    int nlist_ind, double *corr_begin, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  _calc_restricted_point_corr(nlist_ind, ind_list_begin, ind_list_end);
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(
        m_params, m_corr_param_key, *ind_list_begin);
  }
}

/// \brief Calculate select point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_default::_calc_restricted_point_corr(
    int nlist_ind, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  m_params.pre_eval();
  {
    _point_prepare<double>(nlist_ind);
    for (; ind_list_begin < ind_list_end; ind_list_begin++) {
      ParamPack::Val<double>::set(
          m_params, m_corr_param_key, *ind_list_begin,
          (this->*m_flower_func_table_0[nlist_ind][*ind_list_begin])());
    }
  }
  m_params.post_eval();
}

/// \brief Calculate the change in point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_default::_calc_delta_point_corr(
    int nlist_ind, int occ_i, int occ_f, double *corr_begin) const {
  _calc_delta_point_corr(nlist_ind, occ_i, occ_f);
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate the change in point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_default::_calc_delta_point_corr(
    int nlist_ind, int occ_i, int occ_f) const {
  m_params.pre_eval();
  {
    _point_prepare<double>(nlist_ind);
    for (size_type i = 0; i < corr_size(); i++) {
      ParamPack::Val<double>::set(
          m_params, m_corr_param_key, i,
          (this->*m_delta_func_table_0[nlist_ind][i])(occ_i, occ_f));
    }
  }
  m_params.post_eval();
}

/// \brief Calculate the change in select point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_default::_calc_restricted_delta_point_corr(
    int nlist_ind, int occ_i, int occ_f, double *corr_begin,
    size_type const *ind_list_begin, size_type const *ind_list_end) const {
  _calc_restricted_delta_point_corr(nlist_ind, occ_i, occ_f, ind_list_begin,
                                    ind_list_end);
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(
        m_params, m_corr_param_key, *ind_list_begin);
  }
}

/// \brief Calculate the change in select point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_default::_calc_restricted_delta_point_corr(
    int nlist_ind, int occ_i, int occ_f, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  m_params.pre_eval();
  {
    _point_prepare<double>(nlist_ind);
    for (; ind_list_begin < ind_list_end; ind_list_begin++) {
      ParamPack::Val<double>::set(
          m_params, m_corr_param_key, *ind_list_begin,
          (this->*m_delta_func_table_0[nlist_ind][*ind_list_begin])(occ_i,
                                                                    occ_f));
    }
  }
  m_params.post_eval();
}

template <typename Scalar>
void FCC_binary_vacancy_Clexulator_default::_point_prepare(
    int nlist_ind) const {
  switch (nlist_ind) {
    case 0:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 0,
                                    eval_occ_func_0_0(0));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 0,
                                    eval_occ_func_0_1(0));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 1,
                                    eval_occ_func_0_0(1));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 1,
                                    eval_occ_func_0_1(1));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 2,
                                    eval_occ_func_0_0(2));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 2,
                                    eval_occ_func_0_1(2));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 3,
                                    eval_occ_func_0_0(3));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 3,
                                    eval_occ_func_0_1(3));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 4,
                                    eval_occ_func_0_0(4));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 4,
                                    eval_occ_func_0_1(4));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 5,
                                    eval_occ_func_0_0(5));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 5,
                                    eval_occ_func_0_1(5));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 6,
                                    eval_occ_func_0_0(6));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 6,
                                    eval_occ_func_0_1(6));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 7,
                                    eval_occ_func_0_0(7));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 7,
                                    eval_occ_func_0_1(7));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 8,
                                    eval_occ_func_0_0(8));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 8,
                                    eval_occ_func_0_1(8));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 9,
                                    eval_occ_func_0_0(9));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 9,
                                    eval_occ_func_0_1(9));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 10,
                                    eval_occ_func_0_0(10));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 10,
                                    eval_occ_func_0_1(10));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 11,
                                    eval_occ_func_0_0(11));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 11,
                                    eval_occ_func_0_1(11));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 12,
                                    eval_occ_func_0_0(12));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 12,
                                    eval_occ_func_0_1(12));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 13,
                                    eval_occ_func_0_0(13));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 13,
                                    eval_occ_func_0_1(13));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14,
                                    eval_occ_func_0_0(14));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 14,
                                    eval_occ_func_0_1(14));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 15,
                                    eval_occ_func_0_0(15));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 15,
                                    eval_occ_func_0_1(15));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 16,
                                    eval_occ_func_0_0(16));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 16,
                                    eval_occ_func_0_1(16));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 17,
                                    eval_occ_func_0_0(17));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 17,
                                    eval_occ_func_0_1(17));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18,
                                    eval_occ_func_0_0(18));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 18,
                                    eval_occ_func_0_1(18));
      }
      break;
  }
}
template <typename Scalar>
void FCC_binary_vacancy_Clexulator_default::_global_prepare() const {
  if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 0,
                                eval_occ_func_0_0(0));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 0,
                                eval_occ_func_0_1(0));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 1,
                                eval_occ_func_0_0(1));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 1,
                                eval_occ_func_0_1(1));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 2,
                                eval_occ_func_0_0(2));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 2,
                                eval_occ_func_0_1(2));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 3,
                                eval_occ_func_0_0(3));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 3,
                                eval_occ_func_0_1(3));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 4,
                                eval_occ_func_0_0(4));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 4,
                                eval_occ_func_0_1(4));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 5,
                                eval_occ_func_0_0(5));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 5,
                                eval_occ_func_0_1(5));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 6,
                                eval_occ_func_0_0(6));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 6,
                                eval_occ_func_0_1(6));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 7,
                                eval_occ_func_0_0(7));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 7,
                                eval_occ_func_0_1(7));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 8,
                                eval_occ_func_0_0(8));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 8,
                                eval_occ_func_0_1(8));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 9,
                                eval_occ_func_0_0(9));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 9,
                                eval_occ_func_0_1(9));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 10,
                                eval_occ_func_0_0(10));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 10,
                                eval_occ_func_0_1(10));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 11,
                                eval_occ_func_0_0(11));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 11,
                                eval_occ_func_0_1(11));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 12,
                                eval_occ_func_0_0(12));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 12,
                                eval_occ_func_0_1(12));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 13,
                                eval_occ_func_0_0(13));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 13,
                                eval_occ_func_0_1(13));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14,
                                eval_occ_func_0_0(14));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 14,
                                eval_occ_func_0_1(14));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 15,
                                eval_occ_func_0_0(15));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 15,
                                eval_occ_func_0_1(15));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 16,
                                eval_occ_func_0_0(16));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 16,
                                eval_occ_func_0_1(16));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 17,
                                eval_occ_func_0_0(17));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 17,
                                eval_occ_func_0_1(17));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18,
                                eval_occ_func_0_0(18));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 18,
                                eval_occ_func_0_1(18));
  }
}

// Basis functions for empty cluster:
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_0_0() const {
  return 1;
}

/**** Basis functions for orbit 1****
0.0000000 0.0000000 0.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_1_0() const {
  return occ_func_0_0(0);
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_1_1() const {
  return occ_func_0_1(0);
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_1_0_at_0() const {
  return occ_func_0_0(0);
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_1_1_at_0() const {
  return occ_func_0_1(0);
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_1_0_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]);
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_1_1_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]);
}

/**** Basis functions for orbit 2****
0.0000000 0.0000000 0.0000000 A  B  Va

0.0000000 1.0000000 0.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_0() const {
  return (occ_func_0_0(0) * occ_func_0_0(9) +
          occ_func_0_0(10) * occ_func_0_0(0) +
          occ_func_0_0(11) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(12) +
          occ_func_0_0(0) * occ_func_0_0(7) +
          occ_func_0_0(0) * occ_func_0_0(8)) /
         6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_1() const {
  return ((0.707107 * occ_func_0_0(0) * occ_func_0_1(9) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(9)) +
          (0.707107 * occ_func_0_0(10) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(10) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(11) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(11) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(12) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(12)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(7) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(7)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(8) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(8))) /
         6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_2_2() const {
  return (occ_func_0_1(0) * occ_func_0_1(9) +
          occ_func_0_1(10) * occ_func_0_1(0) +
          occ_func_0_1(11) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(12) +
          occ_func_0_1(0) * occ_func_0_1(7) +
          occ_func_0_1(0) * occ_func_0_1(8)) /
         6.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_0_at_0() const {
  return (occ_func_0_0(0) * occ_func_0_0(9) +
          occ_func_0_0(4) * occ_func_0_0(0) +
          occ_func_0_0(10) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(3) +
          occ_func_0_0(11) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(2) +
          occ_func_0_0(0) * occ_func_0_0(12) +
          occ_func_0_0(1) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(7) +
          occ_func_0_0(6) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(8) +
          occ_func_0_0(5) * occ_func_0_0(0)) /
         6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_1_at_0() const {
  return ((0.707107 * occ_func_0_0(0) * occ_func_0_1(9) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(9)) +
          (0.707107 * occ_func_0_0(4) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(4) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(10) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(10) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(3) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(3)) +
          (0.707107 * occ_func_0_0(11) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(11) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(2) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(2)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(12) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(12)) +
          (0.707107 * occ_func_0_0(1) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(1) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(7) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(7)) +
          (0.707107 * occ_func_0_0(6) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(6) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(8) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(8)) +
          (0.707107 * occ_func_0_0(5) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(5) * occ_func_0_0(0))) /
         6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_2_2_at_0() const {
  return (occ_func_0_1(0) * occ_func_0_1(9) +
          occ_func_0_1(4) * occ_func_0_1(0) +
          occ_func_0_1(10) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(3) +
          occ_func_0_1(11) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(2) +
          occ_func_0_1(0) * occ_func_0_1(12) +
          occ_func_0_1(1) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(7) +
          occ_func_0_1(6) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(8) +
          occ_func_0_1(5) * occ_func_0_1(0)) /
         6.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_0_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         (occ_func_0_0(9) + occ_func_0_0(4) + occ_func_0_0(10) +
          occ_func_0_0(3) + occ_func_0_0(11) + occ_func_0_0(2) +
          occ_func_0_0(12) + occ_func_0_0(1) + occ_func_0_0(7) +
          occ_func_0_0(6) + occ_func_0_0(8) + occ_func_0_0(5)) /
         6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_1_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (0.707107 * occ_func_0_1(9) + 0.707107 * occ_func_0_1(4) +
              0.707107 * occ_func_0_1(10) + 0.707107 * occ_func_0_1(3) +
              0.707107 * occ_func_0_1(11) + 0.707107 * occ_func_0_1(2) +
              0.707107 * occ_func_0_1(12) + 0.707107 * occ_func_0_1(1) +
              0.707107 * occ_func_0_1(7) + 0.707107 * occ_func_0_1(6) +
              0.707107 * occ_func_0_1(8) + 0.707107 * occ_func_0_1(5)) /
             6. +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (0.707107 * occ_func_0_0(9) + 0.707107 * occ_func_0_0(4) +
              0.707107 * occ_func_0_0(10) + 0.707107 * occ_func_0_0(3) +
              0.707107 * occ_func_0_0(11) + 0.707107 * occ_func_0_0(2) +
              0.707107 * occ_func_0_0(12) + 0.707107 * occ_func_0_0(1) +
              0.707107 * occ_func_0_0(7) + 0.707107 * occ_func_0_0(6) +
              0.707107 * occ_func_0_0(8) + 0.707107 * occ_func_0_0(5)) /
             6.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_2_2_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         (occ_func_0_1(9) + occ_func_0_1(4) + occ_func_0_1(10) +
          occ_func_0_1(3) + occ_func_0_1(11) + occ_func_0_1(2) +
          occ_func_0_1(12) + occ_func_0_1(1) + occ_func_0_1(7) +
          occ_func_0_1(6) + occ_func_0_1(8) + occ_func_0_1(5)) /
         6.;
}

/**** Basis functions for orbit 3****
0.0000000 0.0000000 0.0000000 A  B  Va

1.0000000 -1.0000000 -1.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_0() const {
  return (occ_func_0_0(0) * occ_func_0_0(16) +
          occ_func_0_0(0) * occ_func_0_0(18) +
          occ_func_0_0(17) * occ_func_0_0(0)) /
         3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_1() const {
  return ((0.707107 * occ_func_0_0(0) * occ_func_0_1(16) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(16)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(18) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(18)) +
          (0.707107 * occ_func_0_0(17) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(17) * occ_func_0_0(0))) /
         3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::eval_bfunc_3_2() const {
  return (occ_func_0_1(0) * occ_func_0_1(16) +
          occ_func_0_1(0) * occ_func_0_1(18) +
          occ_func_0_1(17) * occ_func_0_1(0)) /
         3.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_0_at_0() const {
  return (occ_func_0_0(0) * occ_func_0_0(16) +
          occ_func_0_0(15) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(18) +
          occ_func_0_0(13) * occ_func_0_0(0) +
          occ_func_0_0(17) * occ_func_0_0(0) +
          occ_func_0_0(0) * occ_func_0_0(14)) /
         3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_1_at_0() const {
  return ((0.707107 * occ_func_0_0(0) * occ_func_0_1(16) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(16)) +
          (0.707107 * occ_func_0_0(15) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(15) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(18) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(18)) +
          (0.707107 * occ_func_0_0(13) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(13) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(17) * occ_func_0_1(0) +
           0.707107 * occ_func_0_1(17) * occ_func_0_0(0)) +
          (0.707107 * occ_func_0_0(0) * occ_func_0_1(14) +
           0.707107 * occ_func_0_1(0) * occ_func_0_0(14))) /
         3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_eval_bfunc_3_2_at_0() const {
  return (occ_func_0_1(0) * occ_func_0_1(16) +
          occ_func_0_1(15) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(18) +
          occ_func_0_1(13) * occ_func_0_1(0) +
          occ_func_0_1(17) * occ_func_0_1(0) +
          occ_func_0_1(0) * occ_func_0_1(14)) /
         3.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_0_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
         (occ_func_0_0(16) + occ_func_0_0(15) + occ_func_0_0(18) +
          occ_func_0_0(13) + occ_func_0_0(17) + occ_func_0_0(14)) /
         3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_1_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) *
             (0.707107 * occ_func_0_1(16) + 0.707107 * occ_func_0_1(15) +
              0.707107 * occ_func_0_1(18) + 0.707107 * occ_func_0_1(13) +
              0.707107 * occ_func_0_1(17) + 0.707107 * occ_func_0_1(14)) /
             3. +
         (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
             (0.707107 * occ_func_0_0(16) + 0.707107 * occ_func_0_0(15) +
              0.707107 * occ_func_0_0(18) + 0.707107 * occ_func_0_0(13) +
              0.707107 * occ_func_0_0(17) + 0.707107 * occ_func_0_0(14)) /
             3.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_default::site_deval_bfunc_3_2_at_0(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) *
         (occ_func_0_1(16) + occ_func_0_1(15) + occ_func_0_1(18) +
          occ_func_0_1(13) + occ_func_0_1(17) + occ_func_0_1(14)) /
         3.;
}

}  // namespace clexulator
}  // namespace CASM

extern "C" {
/// \brief Returns a clexulator::BaseClexulator* owning a
/// FCC_binary_vacancy_Clexulator_default
CASM::clexulator::BaseClexulator *make_FCC_binary_vacancy_Clexulator_default() {
  return new CASM::clexulator::FCC_binary_vacancy_Clexulator_default();
}
}
