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
    "method" : "local_max_length",
    "params" : {
      "generating_group" : [ 0, 16, 18, 23, 25, 27, 32, 47 ],
      "orbit_branch_specs" : {
        "0" : {
          "cutoff_radius" : 0.000000000000,
          "max_length" : 0.000000000000
        },
        "1" : {
          "cutoff_radius" : 3.000000000000,
          "max_length" : 0.000000000000
        }
      },
      "phenomenal" : {
        "max_length" : 2.828427124746,
        "min_length" : 2.828427124746,
        "sites" : [
          [ 0, 0, 0, 0 ],
          [ 0, 0, 0, 1 ]
        ]
      }
    }
  }
}

**/

/// \brief Returns a clexulator::BaseClexulator* owning a
/// FCC_binary_vacancy_Clexulator_A_Va_1NN_2
extern "C" CASM::clexulator::BaseClexulator *
make_FCC_binary_vacancy_Clexulator_A_Va_1NN_2();

namespace CASM {
namespace clexulator {

/****** GENERATED CLEXPARAMPACK DEFINITION ******/

typedef BasicClexParamPack ParamPack;

/****** GENERATED CLEXULATOR DEFINITION ******/

class FCC_binary_vacancy_Clexulator_A_Va_1NN_2
    : public clexulator::BaseClexulator {
 public:
  FCC_binary_vacancy_Clexulator_A_Va_1NN_2();

  ~FCC_binary_vacancy_Clexulator_A_Va_1NN_2();

  ClexParamPack const &param_pack() const override { return m_params; }

  ClexParamPack &param_pack() override { return m_params; }

  template <typename Scalar>
  Scalar eval_bfunc_0_0() const;

  template <typename Scalar>
  Scalar eval_bfunc_1_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_1_1() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_1_0_at_5() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_1_1_at_5() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_1_0_at_50() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_1_1_at_50() const;
  template <typename Scalar>
  Scalar site_deval_bfunc_1_0_at_5(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_1_1_at_5(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_1_0_at_50(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_1_1_at_50(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar eval_bfunc_2_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_2_1() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_2_0_at_1() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_1_at_1() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_0_at_14() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_1_at_14() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_0_at_12() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_1_at_12() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_0_at_18() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_2_1_at_18() const;

  template <typename Scalar>
  Scalar site_deval_bfunc_2_0_at_1(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_1_at_1(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_2_0_at_14(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_1_at_14(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_2_0_at_12(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_1_at_12(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_2_0_at_18(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_2_1_at_18(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar eval_bfunc_3_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_3_1() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_3_0_at_3() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_1_at_3() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_0_at_6() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_1_at_6() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_0_at_9() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_1_at_9() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_0_at_11() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_3_1_at_11() const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_0_at_3(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_1_at_3(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_3_0_at_6(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_1_at_6(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_3_0_at_9(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_1_at_9(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_3_0_at_11(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_3_1_at_11(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar eval_bfunc_4_0() const;
  template <typename Scalar>
  Scalar eval_bfunc_4_1() const;

  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_2() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_2() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_26() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_26() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_4() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_4() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_7() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_7() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_31() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_31() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_33() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_33() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_10() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_10() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_0_at_38() const;
  template <typename Scalar>
  Scalar site_eval_bfunc_4_1_at_38() const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_2(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_2(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_26(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_26(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_4(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_4(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_7(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_7(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_31(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_31(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_33(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_33(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_10(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_10(int occ_i, int occ_f) const;

  template <typename Scalar>
  Scalar site_deval_bfunc_4_0_at_38(int occ_i, int occ_f) const;
  template <typename Scalar>
  Scalar site_deval_bfunc_4_1_at_38(int occ_i, int occ_f) const;

 private:
  // ParamPack object, which stores temporary data for calculations
  mutable ParamPack m_params;

  // typedef for method pointers of scalar type double
  typedef double (
      FCC_binary_vacancy_Clexulator_A_Va_1NN_2::*BasisFuncPtr_0)() const;

  // typedef for method pointers
  typedef double (FCC_binary_vacancy_Clexulator_A_Va_1NN_2::*
                      DeltaBasisFuncPtr_0)(int, int) const;

  // array of pointers to member functions for calculating basis functions of
  // scalar type double
  BasisFuncPtr_0 m_orbit_func_table_0[9];

  // array of pointers to member functions for calculating flower functions of
  // scalar type double
  BasisFuncPtr_0 m_flower_func_table_0[51][9];

  // array of pointers to member functions for calculating DELTA flower
  // functions of scalar type double
  DeltaBasisFuncPtr_0 m_delta_func_table_0[51][9];

  // Occupation Function tables for basis sites in asymmetric unit 0:
  //   - basis site 0:
  double m_occ_func_0_0[3];
  double m_occ_func_0_1[3];

  // ClexParamPack allocation for evaluated correlations
  ParamPack::Key m_corr_param_key;
  // ClexParamPack allocation for DoF occ
  ParamPack::Key m_occ_site_func_param_key;

  /// \brief Clone the FCC_binary_vacancy_Clexulator_A_Va_1NN_2
  BaseClexulator *_clone() const override {
    return new FCC_binary_vacancy_Clexulator_A_Va_1NN_2(*this);
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

FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
    FCC_binary_vacancy_Clexulator_A_Va_1NN_2()
    : BaseClexulator(51, 9, 51) {
  m_occ_func_0_0[0] = 0.0000000000, m_occ_func_0_0[1] = 1.0000000000,
  m_occ_func_0_0[2] = 0.0000000000;

  m_occ_func_0_1[0] = 0.0000000000, m_occ_func_0_1[1] = 0.0000000000,
  m_occ_func_0_1[2] = 1.0000000000;

  m_occ_site_func_param_key = m_params.allocate("occ_site_func", 2, 51, true);

  m_corr_param_key = m_params.allocate("corr", corr_size(), 1, false);

  m_orbit_func_table_0[0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_0_0<double>;
  m_orbit_func_table_0[1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_1_0<double>;
  m_orbit_func_table_0[2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_1_1<double>;
  m_orbit_func_table_0[3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_2_0<double>;
  m_orbit_func_table_0[4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_2_1<double>;
  m_orbit_func_table_0[5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_3_0<double>;
  m_orbit_func_table_0[6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_3_1<double>;
  m_orbit_func_table_0[7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_4_0<double>;
  m_orbit_func_table_0[8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_4_1<double>;

  m_flower_func_table_0[0][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[0][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[1][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_1<
          double>;
  m_flower_func_table_0[1][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_1<
          double>;
  m_flower_func_table_0[1][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[1][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[2][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[2][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_2<
          double>;
  m_flower_func_table_0[2][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_2<
          double>;

  m_flower_func_table_0[3][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_3<
          double>;
  m_flower_func_table_0[3][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_3<
          double>;
  m_flower_func_table_0[3][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[3][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[4][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[4][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_4<
          double>;
  m_flower_func_table_0[4][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_4<
          double>;

  m_flower_func_table_0[5][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_0_at_5<
          double>;
  m_flower_func_table_0[5][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_1_at_5<
          double>;
  m_flower_func_table_0[5][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[5][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[6][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_6<
          double>;
  m_flower_func_table_0[6][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_6<
          double>;
  m_flower_func_table_0[6][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[6][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[7][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[7][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_7<
          double>;
  m_flower_func_table_0[7][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_7<
          double>;

  m_flower_func_table_0[8][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[8][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[9][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_9<
          double>;
  m_flower_func_table_0[9][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_9<
          double>;
  m_flower_func_table_0[9][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[9][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[10][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[10][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_10<
          double>;
  m_flower_func_table_0[10][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_10<
          double>;

  m_flower_func_table_0[11][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_11<
          double>;
  m_flower_func_table_0[11][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_11<
          double>;
  m_flower_func_table_0[11][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[11][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[12][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_12<
          double>;
  m_flower_func_table_0[12][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_12<
          double>;
  m_flower_func_table_0[12][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[12][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[13][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[13][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[14][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_14<
          double>;
  m_flower_func_table_0[14][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_14<
          double>;
  m_flower_func_table_0[14][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[14][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[15][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[15][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[16][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[16][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[17][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[17][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[18][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_18<
          double>;
  m_flower_func_table_0[18][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_18<
          double>;
  m_flower_func_table_0[18][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[18][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[19][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[19][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[20][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[20][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[21][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[21][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[22][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[22][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[23][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[23][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[24][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[24][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[25][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[25][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[26][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[26][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_26<
          double>;
  m_flower_func_table_0[26][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_26<
          double>;

  m_flower_func_table_0[27][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[27][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[28][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[28][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[29][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[29][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[30][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[30][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[31][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[31][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_31<
          double>;
  m_flower_func_table_0[31][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_31<
          double>;

  m_flower_func_table_0[32][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[32][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[33][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[33][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_33<
          double>;
  m_flower_func_table_0[33][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_33<
          double>;

  m_flower_func_table_0[34][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[34][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[35][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[35][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[36][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[36][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[37][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[37][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[38][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[38][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_38<
          double>;
  m_flower_func_table_0[38][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_38<
          double>;

  m_flower_func_table_0[39][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[39][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[40][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[40][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[41][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[41][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[42][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[42][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[43][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[43][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[44][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[44][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[45][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[45][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[46][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[46][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[47][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[47][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[48][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[48][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[49][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[49][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_flower_func_table_0[50][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_0_at_50<
          double>;
  m_flower_func_table_0[50][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_1_at_50<
          double>;
  m_flower_func_table_0[50][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_flower_func_table_0[50][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[0][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[0][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[1][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_1<
          double>;
  m_delta_func_table_0[1][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_1<
          double>;
  m_delta_func_table_0[1][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[1][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[2][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[2][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_2<
          double>;
  m_delta_func_table_0[2][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_2<
          double>;

  m_delta_func_table_0[3][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_3<
          double>;
  m_delta_func_table_0[3][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_3<
          double>;
  m_delta_func_table_0[3][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[3][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[4][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[4][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_4<
          double>;
  m_delta_func_table_0[4][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_4<
          double>;

  m_delta_func_table_0[5][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_0_at_5<
          double>;
  m_delta_func_table_0[5][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_1_at_5<
          double>;
  m_delta_func_table_0[5][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[5][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[6][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_6<
          double>;
  m_delta_func_table_0[6][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_6<
          double>;
  m_delta_func_table_0[6][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[6][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[7][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[7][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_7<
          double>;
  m_delta_func_table_0[7][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_7<
          double>;

  m_delta_func_table_0[8][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[8][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[9][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_9<
          double>;
  m_delta_func_table_0[9][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_9<
          double>;
  m_delta_func_table_0[9][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[9][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[10][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[10][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_10<
          double>;
  m_delta_func_table_0[10][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_10<
          double>;

  m_delta_func_table_0[11][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_11<
          double>;
  m_delta_func_table_0[11][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_11<
          double>;
  m_delta_func_table_0[11][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[11][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[12][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_12<
          double>;
  m_delta_func_table_0[12][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_12<
          double>;
  m_delta_func_table_0[12][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[12][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[13][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[13][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[14][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_14<
          double>;
  m_delta_func_table_0[14][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_14<
          double>;
  m_delta_func_table_0[14][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[14][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[15][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[15][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[16][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[16][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[17][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[17][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[18][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_18<
          double>;
  m_delta_func_table_0[18][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_18<
          double>;
  m_delta_func_table_0[18][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[18][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[19][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[19][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[20][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[20][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[21][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[21][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[22][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[22][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[23][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[23][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[24][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[24][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[25][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[25][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[26][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[26][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_26<
          double>;
  m_delta_func_table_0[26][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_26<
          double>;

  m_delta_func_table_0[27][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[27][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[28][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[28][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[29][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[29][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[30][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[30][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[31][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[31][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_31<
          double>;
  m_delta_func_table_0[31][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_31<
          double>;

  m_delta_func_table_0[32][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[32][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[33][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[33][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_33<
          double>;
  m_delta_func_table_0[33][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_33<
          double>;

  m_delta_func_table_0[34][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[34][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[35][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[35][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[36][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[36][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[37][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[37][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[38][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[38][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_38<
          double>;
  m_delta_func_table_0[38][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_38<
          double>;

  m_delta_func_table_0[39][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[39][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[40][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[40][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[41][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[41][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[42][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[42][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[43][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[43][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[44][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[44][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[45][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[45][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[46][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[46][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[47][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[47][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[48][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[48][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[49][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[49][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_delta_func_table_0[50][0] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][1] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_0_at_50<
          double>;
  m_delta_func_table_0[50][2] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_1_at_50<
          double>;
  m_delta_func_table_0[50][3] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][4] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][5] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][6] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][7] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;
  m_delta_func_table_0[50][8] =
      &FCC_binary_vacancy_Clexulator_A_Va_1NN_2::zero_func<double>;

  m_weight_matrix.row(0) << 2, 1, 1;
  m_weight_matrix.row(1) << 1, 2, 1;
  m_weight_matrix.row(2) << 1, 1, 2;

  m_sublat_indices = std::set<int>{0};

  m_n_sublattices = 1;

  m_neighborhood = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, 0, 0),  xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, -1), xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, -1), xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),  xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 1),   xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, 0),   xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, -1),  xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -1),  xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, -2),  xtal::UnitCell(1, 1, -1)};

  m_orbit_neighborhood.resize(corr_size());
  m_orbit_site_neighborhood.resize(corr_size());
  m_orbit_neighborhood[1] = std::set<xtal::UnitCell>{xtal::UnitCell(0, -1, 1),
                                                     xtal::UnitCell(0, 2, -2)};
  m_orbit_neighborhood[2] = m_orbit_neighborhood[1];

  m_orbit_site_neighborhood[1] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, 0, -1, 1), xtal::UnitCellCoord(0, 0, 2, -2)};
  m_orbit_site_neighborhood[2] = m_orbit_site_neighborhood[1];

  m_orbit_neighborhood[3] = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, 0, 0), xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(1, 0, 0), xtal::UnitCell(1, 1, -1)};
  m_orbit_neighborhood[4] = m_orbit_neighborhood[3];

  m_orbit_site_neighborhood[3] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, -1, 0, 0), xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, 1, 0, 0), xtal::UnitCellCoord(0, 1, 1, -1)};
  m_orbit_site_neighborhood[4] = m_orbit_site_neighborhood[3];

  m_orbit_neighborhood[5] = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, 1, 0), xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 1, 0), xtal::UnitCell(1, 0, -1)};
  m_orbit_neighborhood[6] = m_orbit_neighborhood[5];

  m_orbit_site_neighborhood[5] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, -1, 1, 0), xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 1, 0), xtal::UnitCellCoord(0, 1, 0, -1)};
  m_orbit_site_neighborhood[6] = m_orbit_site_neighborhood[5];

  m_orbit_neighborhood[7] = std::set<xtal::UnitCell>{
      xtal::UnitCell(-1, 0, 1), xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(0, -1, 0), xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -2), xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -1, 0), xtal::UnitCell(1, 1, -2)};
  m_orbit_neighborhood[8] = m_orbit_neighborhood[7];

  m_orbit_site_neighborhood[7] = std::set<xtal::UnitCellCoord>{
      xtal::UnitCellCoord(0, -1, 0, 1), xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, 0, -1, 0), xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -2), xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -1, 0), xtal::UnitCellCoord(0, 1, 1, -2)};
  m_orbit_site_neighborhood[8] = m_orbit_site_neighborhood[7];
}

FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
    ~FCC_binary_vacancy_Clexulator_A_Va_1NN_2() {
  // nothing here for now
}

/// \brief Calculate contribution to global correlations from one unit cell
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_global_corr_contribution(
    double *corr_begin) const {
  _calc_global_corr_contribution();
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate contribution to global correlations from one unit cell
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_global_corr_contribution()
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_point_corr(
    int nlist_ind, double *corr_begin) const {
  _calc_point_corr(nlist_ind);
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_point_corr(
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_restricted_point_corr(
    int nlist_ind, double *corr_begin, size_type const *ind_list_begin,
    size_type const *ind_list_end) const {
  _calc_restricted_point_corr(nlist_ind, ind_list_begin, ind_list_end);
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(
        m_params, m_corr_param_key, *ind_list_begin);
  }
}

/// \brief Calculate select point correlations about basis site 'nlist_ind'
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_restricted_point_corr(
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_delta_point_corr(
    int nlist_ind, int occ_i, int occ_f, double *corr_begin) const {
  _calc_delta_point_corr(nlist_ind, occ_i, occ_f);
  for (size_type i = 0; i < corr_size(); i++) {
    *(corr_begin + i) =
        ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
  }
}

/// \brief Calculate the change in point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_calc_delta_point_corr(
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
    _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f,
                                      double *corr_begin,
                                      size_type const *ind_list_begin,
                                      size_type const *ind_list_end) const {
  _calc_restricted_delta_point_corr(nlist_ind, occ_i, occ_f, ind_list_begin,
                                    ind_list_end);
  for (; ind_list_begin < ind_list_end; ind_list_begin++) {
    *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(
        m_params, m_corr_param_key, *ind_list_begin);
  }
}

/// \brief Calculate the change in select point correlations due to changing an
/// occupant
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::
    _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f,
                                      size_type const *ind_list_begin,
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
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_point_prepare(
    int nlist_ind) const {
  switch (nlist_ind) {
    case 1:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 1,
                                    eval_occ_func_0_0(1));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 1,
                                    eval_occ_func_0_1(1));
      }
      break;
    case 2:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 2,
                                    eval_occ_func_0_0(2));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 2,
                                    eval_occ_func_0_1(2));
      }
      break;
    case 14:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14,
                                    eval_occ_func_0_0(14));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 14,
                                    eval_occ_func_0_1(14));
      }
      break;
    case 3:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 3,
                                    eval_occ_func_0_0(3));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 3,
                                    eval_occ_func_0_1(3));
      }
      break;
    case 26:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 26,
                                    eval_occ_func_0_0(26));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 26,
                                    eval_occ_func_0_1(26));
      }
      break;
    case 4:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 4,
                                    eval_occ_func_0_0(4));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 4,
                                    eval_occ_func_0_1(4));
      }
      break;
    case 5:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 5,
                                    eval_occ_func_0_0(5));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 5,
                                    eval_occ_func_0_1(5));
      }
      break;
    case 6:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 6,
                                    eval_occ_func_0_0(6));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 6,
                                    eval_occ_func_0_1(6));
      }
      break;
    case 7:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 7,
                                    eval_occ_func_0_0(7));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 7,
                                    eval_occ_func_0_1(7));
      }
      break;
    case 31:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 31,
                                    eval_occ_func_0_0(31));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 31,
                                    eval_occ_func_0_1(31));
      }
      break;
    case 9:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 9,
                                    eval_occ_func_0_0(9));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 9,
                                    eval_occ_func_0_1(9));
      }
      break;
    case 50:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 50,
                                    eval_occ_func_0_0(50));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 50,
                                    eval_occ_func_0_1(50));
      }
      break;
    case 33:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 33,
                                    eval_occ_func_0_0(33));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 33,
                                    eval_occ_func_0_1(33));
      }
      break;
    case 10:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 10,
                                    eval_occ_func_0_0(10));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 10,
                                    eval_occ_func_0_1(10));
      }
      break;
    case 11:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 11,
                                    eval_occ_func_0_0(11));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 11,
                                    eval_occ_func_0_1(11));
      }
      break;
    case 12:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 12,
                                    eval_occ_func_0_0(12));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 12,
                                    eval_occ_func_0_1(12));
      }
      break;
    case 38:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 38,
                                    eval_occ_func_0_0(38));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 38,
                                    eval_occ_func_0_1(38));
      }
      break;
    case 18:
      if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18,
                                    eval_occ_func_0_0(18));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 18,
                                    eval_occ_func_0_1(18));
      }
      break;
  }
}
template <typename Scalar>
void FCC_binary_vacancy_Clexulator_A_Va_1NN_2::_global_prepare() const {
  if (m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
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
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14,
                                eval_occ_func_0_0(14));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 14,
                                eval_occ_func_0_1(14));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18,
                                eval_occ_func_0_0(18));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 18,
                                eval_occ_func_0_1(18));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 26,
                                eval_occ_func_0_0(26));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 26,
                                eval_occ_func_0_1(26));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 31,
                                eval_occ_func_0_0(31));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 31,
                                eval_occ_func_0_1(31));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 33,
                                eval_occ_func_0_0(33));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 33,
                                eval_occ_func_0_1(33));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 38,
                                eval_occ_func_0_0(38));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 38,
                                eval_occ_func_0_1(38));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 50,
                                eval_occ_func_0_0(50));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 1, 50,
                                eval_occ_func_0_1(50));
  }
}

// Basis functions for empty cluster:
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_0_0() const {
  return 1;
}

/**** Basis functions for orbit 1****
0.0000000 -1.0000000 1.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_1_0() const {
  return (occ_func_0_0(5) + occ_func_0_0(50)) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_1_1() const {
  return (occ_func_0_1(5) + occ_func_0_1(50)) / 2.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_0_at_5()
    const {
  return (occ_func_0_0(5)) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_1_at_5()
    const {
  return (occ_func_0_1(5)) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_0_at_50()
    const {
  return (occ_func_0_0(50)) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_1_1_at_50()
    const {
  return (occ_func_0_1(50)) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_0_at_5(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_1_at_5(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 2.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_0_at_50(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 2.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_1_1_at_50(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 2.;
}

/**** Basis functions for orbit 2****
-1.0000000 0.0000000 0.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_2_0() const {
  return (occ_func_0_0(1) + occ_func_0_0(14) + occ_func_0_0(12) +
          occ_func_0_0(18)) /
         4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_2_1() const {
  return (occ_func_0_1(1) + occ_func_0_1(14) + occ_func_0_1(12) +
          occ_func_0_1(18)) /
         4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_1()
    const {
  return (occ_func_0_0(1)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_1()
    const {
  return (occ_func_0_1(1)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_14()
    const {
  return (occ_func_0_0(14)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_14()
    const {
  return (occ_func_0_1(14)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_12()
    const {
  return (occ_func_0_0(12)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_12()
    const {
  return (occ_func_0_1(12)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_0_at_18()
    const {
  return (occ_func_0_0(18)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_2_1_at_18()
    const {
  return (occ_func_0_1(18)) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_1(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_1(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_14(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_14(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_12(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_12(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_0_at_18(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_2_1_at_18(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

/**** Basis functions for orbit 3****
-1.0000000 1.0000000 0.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_3_0() const {
  return (occ_func_0_0(3) + occ_func_0_0(6) + occ_func_0_0(11) +
          occ_func_0_0(9)) /
         4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_3_1() const {
  return (occ_func_0_1(3) + occ_func_0_1(6) + occ_func_0_1(11) +
          occ_func_0_1(9)) /
         4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_3()
    const {
  return (occ_func_0_0(3)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_3()
    const {
  return (occ_func_0_1(3)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_6()
    const {
  return (occ_func_0_0(6)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_6()
    const {
  return (occ_func_0_1(6)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_9()
    const {
  return (occ_func_0_0(9)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_9()
    const {
  return (occ_func_0_1(9)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_0_at_11()
    const {
  return (occ_func_0_0(11)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_3_1_at_11()
    const {
  return (occ_func_0_1(11)) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_3(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_3(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_6(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_6(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_9(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_9(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_0_at_11(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 4.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_3_1_at_11(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 4.;
}

/**** Basis functions for orbit 4****
-1.0000000 0.0000000 1.0000000 A  B  Va

****/
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_4_0() const {
  return (occ_func_0_0(2) + occ_func_0_0(31) + occ_func_0_0(10) +
          occ_func_0_0(33) + occ_func_0_0(7) + occ_func_0_0(26) +
          occ_func_0_0(4) + occ_func_0_0(38)) /
         8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::eval_bfunc_4_1() const {
  return (occ_func_0_1(2) + occ_func_0_1(31) + occ_func_0_1(10) +
          occ_func_0_1(33) + occ_func_0_1(7) + occ_func_0_1(26) +
          occ_func_0_1(4) + occ_func_0_1(38)) /
         8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_2()
    const {
  return (occ_func_0_0(2)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_2()
    const {
  return (occ_func_0_1(2)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_26()
    const {
  return (occ_func_0_0(26)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_26()
    const {
  return (occ_func_0_1(26)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_4()
    const {
  return (occ_func_0_0(4)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_4()
    const {
  return (occ_func_0_1(4)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_7()
    const {
  return (occ_func_0_0(7)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_7()
    const {
  return (occ_func_0_1(7)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_31()
    const {
  return (occ_func_0_0(31)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_31()
    const {
  return (occ_func_0_1(31)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_33()
    const {
  return (occ_func_0_0(33)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_33()
    const {
  return (occ_func_0_1(33)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_10()
    const {
  return (occ_func_0_0(10)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_10()
    const {
  return (occ_func_0_1(10)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_0_at_38()
    const {
  return (occ_func_0_0(38)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_eval_bfunc_4_1_at_38()
    const {
  return (occ_func_0_1(38)) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_2(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_2(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_26(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_26(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_4(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_4(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_7(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_7(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_31(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_31(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_33(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_33(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_10(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_10(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_0_at_38(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (1) / 8.;
}
template <typename Scalar>
Scalar FCC_binary_vacancy_Clexulator_A_Va_1NN_2::site_deval_bfunc_4_1_at_38(
    int occ_i, int occ_f) const {
  return (m_occ_func_0_1[occ_f] - m_occ_func_0_1[occ_i]) * (1) / 8.;
}

}  // namespace clexulator
}  // namespace CASM

extern "C" {
/// \brief Returns a clexulator::BaseClexulator* owning a
/// FCC_binary_vacancy_Clexulator_A_Va_1NN_2
CASM::clexulator::BaseClexulator *
make_FCC_binary_vacancy_Clexulator_A_Va_1NN_2() {
  return new CASM::clexulator::FCC_binary_vacancy_Clexulator_A_Va_1NN_2();
}
}
