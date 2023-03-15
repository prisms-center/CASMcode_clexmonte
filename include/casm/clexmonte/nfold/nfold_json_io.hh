#ifndef CASM_clexmonte_nfold_json_io
#define CASM_clexmonte_nfold_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/nfold/nfold_impl.hh"

namespace CASM {
namespace clexmonte {
namespace nfold {

template <typename EngineType>
void parse(InputParser<Nfold<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  parser.value =
      std::make_unique<Nfold<EngineType>>(system, random_number_engine);
}

}  // namespace nfold
}  // namespace clexmonte
}  // namespace CASM

#endif
