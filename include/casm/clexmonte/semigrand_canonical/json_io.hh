#ifndef CASM_clexmonte_semigrand_canonical_json_io
#define CASM_clexmonte_semigrand_canonical_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/semigrand_canonical/calculator.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexmonte/state/io/json/parse_conditions_impl.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

template <typename EngineType>
void parse(InputParser<SemiGrandCanonical<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  parser.value = std::make_unique<SemiGrandCanonical<EngineType>>(system);
}

inline void parse(InputParser<SemiGrandCanonicalConditions> &parser,
                  std::shared_ptr<system_type> system, bool is_increment) {
  double temperature_is_zero_tol = 1e-10; /*TODO*/
  parser.value = std::make_unique<SemiGrandCanonicalConditions>(
      get_composition_converter(*system), temperature_is_zero_tol);

  parse_temperature(parser);

  double param_chem_pot_tol = CASM::TOL; /*TODO*/
  parse_param_chem_pot(parser, get_composition_converter(*system));
}

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
