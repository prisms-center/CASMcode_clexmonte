#ifndef CASM_clexmonte_canonical_json_io
#define CASM_clexmonte_canonical_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/canonical/canonical_impl.hh"

namespace CASM {
namespace clexmonte {
namespace canonical {

template <typename EngineType>
void parse(InputParser<Canonical<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  // currently no options
  parser.value = std::make_unique<Canonical<EngineType>>(system);
}

}  // namespace canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
