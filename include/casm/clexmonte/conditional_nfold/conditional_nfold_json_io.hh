#ifndef CASM_clexmonte_conditional_nfold_json_io
#define CASM_clexmonte_conditional_nfold_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/conditional_nfold/conditional_nfold.hh"

namespace CASM {
namespace clexmonte {
namespace conditional_nfold {

template <typename EngineType>
void parse(InputParser<ConditionalNfold<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  double nfold_begin_acceptance_rate = 0.0;
  parser.optional(nfold_begin_acceptance_rate, "nfold_begin_acceptance_rate");

  double nfold_end_acceptance_rate = 1.0;
  parser.optional(nfold_end_acceptance_rate, "nfold_end_acceptance_rate");

  if (!parser.valid()) {
    return;
  }

  parser.value = std::make_unique<ConditionalNfold<EngineType>>(
      system, random_number_engine, nfold_begin_acceptance_rate,
      nfold_end_acceptance_rate);
}

}  // namespace conditional_nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
