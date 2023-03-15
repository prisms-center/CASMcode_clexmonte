#ifndef CASM_clexmonte_kinetic_json_io
#define CASM_clexmonte_kinetic_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/kinetic/kinetic.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

template <typename EngineType>
void parse(InputParser<Kinetic<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  // currently no options
  parser.value =
      std::make_unique<Kinetic<EngineType>>(system, random_number_engine);
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
