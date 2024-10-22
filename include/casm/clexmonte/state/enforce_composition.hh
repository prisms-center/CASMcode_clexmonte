#ifndef CASM_clexmonte_state_enforce_composition
#define CASM_clexmonte_state_enforce_composition

#include <vector>

#include "casm/composition/CompositionCalculator.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/global/eigen.hh"
#include "casm/misc/Validator.hh"
#include "casm/misc/algorithm.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {

namespace composition {
class CompositionCalculator;
}

namespace monte {
class OccSwap;
class OccLocation;
}  // namespace monte

namespace clexmonte {

/// \brief Apply grand canonical swaps to enforce target composition
template <typename GeneratorType>
void enforce_composition(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_mol_composition,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<monte::OccSwap> const &semigrand_canonical_swaps,
    monte::OccLocation &occ_location, GeneratorType &random_number_generator);

/// \brief Enforce composition conditions consistency
void enforce_composition_consistency(
    state_type &state,
    composition::CompositionConverter const &composition_converter,
    double mol_composition_tol);

/// \brief Validate composition conditions consistency
Validator validate_composition_consistency(
    state_type &state,
    composition::CompositionConverter const &composition_converter,
    double mol_composition_tol);

// --- Implementation ---

namespace enforce_composition_impl {

template <typename GeneratorType>
std::vector<monte::OccSwap>::const_iterator find_semigrand_canonical_swap(
    Eigen::VectorXi &occupation, Eigen::VectorXd &current_mol_composition,
    Eigen::VectorXd const &target_mol_composition,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<Index> const &species_to_component_index_converter,
    GeneratorType &random_number_generator,
    monte::OccLocation const &occ_location,
    std::vector<monte::OccSwap>::const_iterator begin,
    std::vector<monte::OccSwap>::const_iterator end) {
  auto const &index_converter = species_to_component_index_converter;

  double original_dist =
      (current_mol_composition - target_mol_composition).norm();
  double best_dist = original_dist;

  double volume = occupation.size() / composition_calculator.n_sublat();
  double dn = 1. / volume;
  double tol = dn * 1e-3;

  // store <distance_to_target_mol_composition>:{swap_iterator, number of
  // swaps}
  typedef std::vector<monte::OccSwap>::const_iterator iterator_type;
  typedef std::tuple<iterator_type, Index, Eigen::VectorXd> value_type;
  std::vector<value_type> choices;

  // check each possible swap for how close the composition is afterwards
  for (auto it = begin; it != end; ++it) {
    if (occ_location.cand_size(it->cand_a)) {
      Eigen::VectorXd tmol_composition = current_mol_composition;
      tmol_composition[index_converter[it->cand_a.species_index]] -= dn;
      tmol_composition[index_converter[it->cand_b.species_index]] += dn;
      double dist = (tmol_composition - target_mol_composition).norm();

      // if no clear improvement, skip
      if (dist > original_dist - tol) {
        continue;

        // if clear improvement, new best
      } else if (dist < best_dist - tol) {
        choices.clear();
        choices.push_back(
            {it, occ_location.cand_size(it->cand_a), tmol_composition});
        best_dist = dist;

        // if tied with existing improvement, add as a choice
      } else if (dist < best_dist + tol) {
        choices.push_back(
            {it, occ_location.cand_size(it->cand_a), tmol_composition});
      }
    }
  }
  if (!choices.size()) {
    return end;
  }

  // break ties randomly, weighted by number of candidates
  double sum = 0.0;
  for (const auto &val : choices) {
    sum += std::get<1>(val);
  }

  double rand = random_number_generator.random_real(sum);
  sum = 0.0;
  for (const auto &val : choices) {
    sum += std::get<1>(val);
    if (rand < sum) {
      current_mol_composition = std::get<2>(val);
      return std::get<0>(val);
    }
  }
  throw std::runtime_error(
      "Error in CASM::clexmonte::find_semigrand_canonical_swap, failed "
      "enforcing "
      "composition");
};

inline std::vector<Index> make_species_to_component_index_converter(
    composition::CompositionCalculator const &composition_calculator,
    monte::Conversions const &convert) {
  std::string msg =
      "Error in CASM::clexmonte::enforce_composition: inconsistency between "
      "composition_calculator and index converter";

  auto const &components = composition_calculator.components();
  Index species_size = convert.species_size();
  if (species_size != components.size()) {
    throw std::runtime_error(msg);
  }

  std::vector<Index> species_to_component_index_converter(species_size);
  auto begin = components.begin();
  auto end = components.end();
  for (Index i_species = 0; i_species < species_size; ++i_species) {
    auto it = std::find(begin, end, convert.species_name(i_species));
    if (it == end) {
      throw std::runtime_error(msg);
    }
    species_to_component_index_converter[i_species] = std::distance(begin, it);
  }
  return species_to_component_index_converter;
}

}  // namespace enforce_composition_impl

/// \brief Apply grand canonical swaps to enforce target composition
///
/// \param occupation Current occupation
/// \param occ_location Current occupant location tracking. Must already be
/// initialized.
///
/// Method:
/// - Find which of the provided grand canonical swap types transforms
///   the composition most closely to the target composition
/// - If no swap can improve the composition, return
/// - Propose and apply an event consistent with the found swap type
/// - Repeat
template <typename GeneratorType>
void enforce_composition(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_mol_composition,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<monte::OccSwap> const &semigrand_canonical_swaps,
    monte::OccLocation &occ_location, GeneratorType &random_number_generator) {
  monte::Conversions const &convert = occ_location.convert();

  // no guarantee convert species_index corresponds to mol_composition index
  std::vector<Index> species_to_component_index_converter =
      enforce_composition_impl::make_species_to_component_index_converter(
          composition_calculator, convert);

  auto begin = semigrand_canonical_swaps.begin();
  auto end = semigrand_canonical_swaps.end();
  monte::OccEvent event;
  Eigen::VectorXd current_mol_composition =
      composition_calculator.mean_num_each_component(occupation);

  while (true) {
    auto it = enforce_composition_impl::find_semigrand_canonical_swap(
        occupation, current_mol_composition, target_mol_composition,
        composition_calculator, species_to_component_index_converter,
        random_number_generator, occ_location, begin, end);

    if (it == end) {
      break;
    }

    /// propose event of chosen candidate type and apply swap
    monte::propose_semigrand_canonical_event_from_swap(event, occ_location, *it,
                                                       random_number_generator);
    occ_location.apply(event, occupation);
  }
}

/// \brief Enforce composition conditions consistency
///
/// - If both present and not consistent, set param_composition to be
///   consistent with mol_composition and print warning
/// - If only one set, set the other to be consistent
inline void enforce_composition_consistency(
    state_type &state,
    composition::CompositionConverter const &composition_converter,
    double mol_composition_tol) {
  monte::ValueMap const &conditions = state.conditions;
  if (conditions.vector_values.count("mol_composition") &&
      conditions.vector_values.count("param_composition")) {
    Eigen::VectorXd mol_composition =
        conditions.vector_values.at("mol_composition");
    Eigen::VectorXd param_composition =
        conditions.vector_values.at("param_composition");
    Eigen::VectorXd equiv_mol_composition =
        composition_converter.mol_composition(param_composition);
    if (!CASM::almost_equal(mol_composition, equiv_mol_composition,
                            mol_composition_tol)) {
      auto &log = CASM::log();
      log.warning<Log::quiet>("Composition conditions mismatch");
      log.indent() << "mol_composition conditions are not consistent with "
                      "param_composition conditions!"
                   << std::endl;
      log.indent() << "mol_composition: " << mol_composition.transpose()
                   << std::endl;
      log.indent() << "param_composition: " << param_composition.transpose()
                   << std::endl;
      log.indent() << "equivalent mol_composition: "
                   << equiv_mol_composition.transpose() << std::endl;
      log.indent() << "Will proceed using mol_composition" << std::endl
                   << std::endl;
    }
  } else if (conditions.vector_values.count("mol_composition")) {
    Eigen::VectorXd mol_composition =
        conditions.vector_values.at("mol_composition");
    Eigen::VectorXd equiv_param_composition =
        composition_converter.param_composition(mol_composition);
    state.conditions.vector_values["param_composition"] =
        equiv_param_composition;
  } else if (conditions.vector_values.count("param_composition")) {
    Eigen::VectorXd param_composition =
        conditions.vector_values.at("param_composition");
    Eigen::VectorXd equiv_mol_composition =
        composition_converter.mol_composition(param_composition);
    state.conditions.vector_values["mol_composition"] = equiv_mol_composition;
  } else {
    throw std::runtime_error(
        "Error in enforce_composition_consistency: Conditions must include "
        "`mol_composition` or `param_composition` (or both).");
  }
}

/// \brief Validate composition conditions consistency
inline Validator validate_composition_consistency(
    state_type &state,
    composition::CompositionConverter const &composition_converter,
    double mol_composition_tol) {
  Validator validator;
  monte::ValueMap const &conditions = state.conditions;
  if (conditions.vector_values.count("mol_composition") &&
      conditions.vector_values.count("param_composition")) {
    Eigen::VectorXd mol_composition =
        conditions.vector_values.at("mol_composition");
    Eigen::VectorXd param_composition =
        conditions.vector_values.at("param_composition");
    Eigen::VectorXd equiv_mol_composition =
        composition_converter.mol_composition(param_composition);
    if (!CASM::almost_equal(mol_composition, equiv_mol_composition,
                            mol_composition_tol)) {
      auto &log = CASM::log();
      std::stringstream msg;
      msg << "mol_composition conditions are not consistent with "
             "param_composition conditions."
          << " mol_composition: " << mol_composition.transpose() << ";"
          << " param_composition: " << param_composition.transpose() << ";"
          << " equivalent mol_composition: "
          << equiv_mol_composition.transpose();

      validator.error.insert(msg.str());
    }
    return validator;
  } else if (conditions.vector_values.count("mol_composition")) {
    return validator;
  } else if (conditions.vector_values.count("param_composition")) {
    return validator;
  } else {
    std::stringstream msg;
    msg << "Neither `mol_composition` nor `param_composition` is included.";
    validator.error.insert(msg.str());
    return validator;
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
