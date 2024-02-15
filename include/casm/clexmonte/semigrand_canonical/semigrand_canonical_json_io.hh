#ifndef CASM_clexmonte_semigrand_canonical_json_io
#define CASM_clexmonte_semigrand_canonical_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/semigrand_canonical/semigrand_canonical.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

template <typename EngineType>
void parse(InputParser<SemiGrandCanonical<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  parser.value = std::make_unique<SemiGrandCanonical<EngineType>>(
      system, random_number_engine);
}

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
