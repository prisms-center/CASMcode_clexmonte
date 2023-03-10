#include "casm/clexmonte/kinetic/kinetic_impl.hh"
#include "casm/clexmonte/state/make_conditions.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Construct a list of atom names corresponding to OccLocation atoms
///
/// Notes:
/// - If atoms are conserved, then the order of this list will remain unchanged
///   during the course of a calculation
/// - Values are set to -1 if atom is no longer in supercell
std::vector<Index> make_atom_name_index_list(
    monte::OccLocation const &occ_location,
    occ_events::OccSystem const &occ_system) {
  // sanity check:
  monte::Conversions const &convert = occ_location.convert();
  if (convert.species_size() != occ_system.orientation_name_list.size()) {
    throw std::runtime_error(
        "Error in CASM::clexmonte::kinetic::make_snapshot_for_conserved_atoms: "
        "mismatch between monte::Conversions and occ_events::OccSystem.");
  }

  // collect atom name indices
  std::vector<Index> atom_name_index_list(occ_location.atom_size(), -1);
  for (Index i = 0; i < occ_location.mol_size(); ++i) {
    monte::Mol const &mol = occ_location.mol(i);
    Index b = convert.l_to_b(mol.l);
    Index occupant_index =
        occ_system.orientation_to_occupant_index[b][mol.species_index];
    Index atom_position_index = 0;
    for (Index atom_id : mol.component) {
      Index atom_name_index =
          occ_system.atom_position_to_name_index[b][occupant_index]
                                                [atom_position_index];
      atom_name_index_list.at(atom_id) = atom_name_index;
      ++atom_position_index;
    }
  }
  return atom_name_index_list;
}

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations
///
/// \param temperature The temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert between species per unit cell
///     (mol_composition) and parameteric composition (param_composition).
/// \param comp A map of component names (for species per unit cell composition)
/// or
///     axes names (for parametric composition) to value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "mol_composition".
///
/// Example: Specifying "mol_composition"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 2.0},            // composition values (#/unit cell)
///     {"O", 1./6.},
///     {"Va", 5./6.}});
/// \endcode
///
/// Example: Specifying "param_composition"
/// \code
/// ValueMap conditions = canonical::make_conditions(
///    300.0,                   // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 1./6.}});         // composition values (param_composition)
/// \endcode
///
monte::ValueMap make_conditions(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["mol_composition"] =
      make_mol_composition(composition_converter, comp);
  return conditions;
}

/// \brief Helper for making a conditions ValueMap for canonical Monte
///     Carlo calculations, interpreted as an increment
///
/// \param temperature The change in temperature
/// \param composition_converter composition::CompositionConverter, used to
///     validate input and convert between species per unit cell
///     (mol_composition) and parameteric composition (param_composition).
/// \param comp A map of component names (for change in mol per unit cell
///     composition) or axes names (for change in parametric composition) to
///     value.
///
/// \returns ValueMap which contains scalar "temperature" and vector
///     "mol_composition" (increment).
///
/// Example: Specifying "mol_composition" increment
/// \code
/// ValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"Zr", 0.0},            // composition values (#/unit cell)
///     {"O", 0.01},
///     {"Va", -0.01}});
/// \endcode
///
/// Example: Specifying "param_composition" increment
/// \code
/// ValueMap conditions_increment = canonical::make_conditions_increment(
///    10.0,                    // temperature (K)
///    composition_converter,   // composition converter
///    {{"a", 0.02}});          // composition values (param_composition)
/// \endcode
///
monte::ValueMap make_conditions_increment(
    double temperature,
    composition::CompositionConverter const &composition_converter,
    std::map<std::string, double> comp) {
  monte::ValueMap conditions;
  conditions.scalar_values["temperature"] = temperature;
  conditions.vector_values["mol_composition"] =
      make_mol_composition_increment(composition_converter, comp);
  return conditions;
}

template struct Kinetic<std::mt19937_64>;

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM
