#ifndef CASM_clexmonte_nfold
#define CASM_clexmonte_nfold

#include "casm/clexmonte/nfold/nfold_events.hh"
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical.hh"
#include "casm/monte/methods/nfold.hh"

namespace CASM {
namespace clexmonte {
namespace nfold {

using semi_grand_canonical::make_conditions;
using semi_grand_canonical::make_conditions_increment;

/// \brief Implements semi-grand canonical Monte Carlo calculations
template <typename EngineType>
struct Nfold : public semi_grand_canonical::SemiGrandCanonical<EngineType> {
  explicit Nfold(std::shared_ptr<system_type> _system,
                 std::shared_ptr<EngineType> _random_number_engine =
                     std::shared_ptr<EngineType>());

  /// Data for N-fold way implementation
  std::shared_ptr<NfoldEventData> event_data;

  /// Data for sampling functions
  monte::NfoldData<config_type> nfold_data;

  /// \brief Perform a single run, evolving current state
  void run(state_type &state, monte::OccLocation &occ_location,
           run_manager_type &run_manager);

  typedef semi_grand_canonical::SemiGrandCanonical<EngineType> Base;
  using Base::standard_analysis_functions;
  using Base::standard_modifying_functions;
  using Base::standard_sampling_functions;
};

/// \brief Explicitly instantiated Nfold calculator
typedef Nfold<std::mt19937_64> Nfold_mt19937_64;

}  // namespace nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
