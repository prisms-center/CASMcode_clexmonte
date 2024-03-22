#ifndef CASM_clexmonte_nfold
#define CASM_clexmonte_nfold

#include "casm/clexmonte/nfold/nfold_events.hh"
#include "casm/clexmonte/semigrand_canonical/calculator.hh"
#include "casm/monte/methods/nfold.hh"

namespace CASM {
namespace clexmonte {
namespace nfold {

using semigrand_canonical::make_conditions;
using semigrand_canonical::make_conditions_increment;

/// \brief Implements semi-grand canonical Monte Carlo calculations
template <typename EngineType>
struct Nfold : public semigrand_canonical::SemiGrandCanonical<EngineType> {
  typedef EngineType engine_type;

  explicit Nfold(std::shared_ptr<system_type> _system);

  /// Method allows time-based sampling
  bool time_sampling_allowed = true;

  /// Data for N-fold way implementation
  std::shared_ptr<NfoldEventData> event_data;

  /// Data for sampling functions
  monte::NfoldData<config_type, statistics_type, engine_type> nfold_data;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type<EngineType> &run_manager);

  typedef semigrand_canonical::SemiGrandCanonical<EngineType> Base;
  using Base::standard_analysis_functions;
  using Base::standard_json_sampling_functions;
  using Base::standard_modifying_functions;
  using Base::standard_sampling_functions;
};

/// \brief Explicitly instantiated Nfold calculator
typedef Nfold<std::mt19937_64> Nfold_mt19937_64;

}  // namespace nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
