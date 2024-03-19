#ifndef CASM_clexmonte_semigrand_canonical_json_io
#define CASM_clexmonte_semigrand_canonical_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/semigrand_canonical/semigrand_canonical.hh"
#include "casm/clexmonte/semigrand_canonical/semigrand_canonical_conditions.hh"

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

inline void parse(InputParser<SemiGrandCanonicalConditions> &parser,
                  std::shared_ptr<system_type> system, bool is_increment) {
  parser.value = std::make_unique<SemiGrandCanonicalConditions>(
      get_composition_converter(*system), );

  parse_temperature(parser);

  double param_chem_pot_tol = CASM::TOL; /*TODO*/
  parse_param_chem_pot(parser, get_composition_converter(*system), is_increment,
                       param_chem_pot_tol);
}

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
