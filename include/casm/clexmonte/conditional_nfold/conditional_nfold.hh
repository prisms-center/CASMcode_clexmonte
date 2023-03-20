#ifndef CASM_clexmonte_conditional_nfold
#define CASM_clexmonte_conditional_nfold

#include "casm/clexmonte/nfold/nfold.hh"
#include "casm/clexmonte/nfold/nfold_events.hh"

namespace CASM {
namespace clexmonte {
namespace conditional_nfold {

using semi_grand_canonical::make_conditions;
using semi_grand_canonical::make_conditions_increment;

/// \brief Implements semi-grand canonical Monte Carlo calculations
template <typename EngineType>
struct ConditionalNfold : public nfold::Nfold<EngineType> {
  explicit ConditionalNfold(std::shared_ptr<system_type> _system,
                            std::shared_ptr<EngineType> _random_number_engine =
                                std::shared_ptr<EngineType>(),
                            double _nfold_begin_acceptance_rate = 0.0,
                            double _nfold_end_acceptance_rate = 1.0);

  double nfold_begin_acceptance_rate;
  double nfold_end_acceptance_rate;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<EngineType> &run_manager);

  typedef nfold::Nfold<EngineType> Base;
  using Base::standard_analysis_functions;
  using Base::standard_modifying_functions;
  using Base::standard_sampling_functions;
};

/// \brief Explicitly instantiated ConditionalNfold calculator
typedef ConditionalNfold<std::mt19937_64> ConditionalNfold_mt19937_64;

}  // namespace conditional_nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
