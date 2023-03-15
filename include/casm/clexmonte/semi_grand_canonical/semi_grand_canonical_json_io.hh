#ifndef CASM_clexmonte_semi_grand_canonical_json_io
#define CASM_clexmonte_semi_grand_canonical_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/semi_grand_canonical/semi_grand_canonical.hh"

namespace CASM {
namespace clexmonte {
namespace semi_grand_canonical {

template <typename EngineType>
void parse(InputParser<SemiGrandCanonical<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  parser.value = std::make_unique<SemiGrandCanonical<EngineType>>(
      system, random_number_engine);
}

}  // namespace semi_grand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
