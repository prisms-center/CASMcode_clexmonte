#ifndef CASM_clexmonte_semigrand_canonical_conditions
#define CASM_clexmonte_semigrand_canonical_conditions

#include "casm/clexmonte/state/Conditions.hh"
#include "casm/monte/ValueMap.hh"

namespace CASM {

namespace composition {
class CompositionConverter;
}

namespace clexmonte {
namespace semigrand_canonical {

/// \brief Holds conditions in form preferable to monte::ValueMap for
/// calculation
///
/// Notes:
/// - Can also be used to specify a conditions increment when specifying a path
///   in parameter space
/// - Can be converted to/from a monte::ValueMap which is more convenient for
///   incrementing, etc.
struct SemiGrandCanonicalConditions : public TemperatureConditionsMixin,
                                      public ParamChemPotConditionsMixin {
  /// \brief Construct default SemiGrandCanonicalConditions
  SemiGrandCanonicalConditions(
      composition::CompositionConverter const &_composition_converter,
      double _temperature_is_zero_tol = 1e-10)
      : TemperatureConditionsMixin(_temperature_is_zero_tol),
        ParamChemPotConditionsMixin(_composition_converter) {}

  void set_all(monte::ValueMap const &map, bool is_increment) {
    this->set_temperature(map);
    this->set_param_chem_pot(map);
  }

  monte::ValueMap to_value_map(bool is_increment) const {
    monte::ValueMap map;
    this->put_temperature(map);
    this->put_param_chem_pot(map);
    return map;
  }
};

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
