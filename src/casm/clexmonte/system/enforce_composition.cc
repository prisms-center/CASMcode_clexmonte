#include "casm/clexmonte/system/enforce_composition.hh"

#include "casm/composition/CompositionCalculator.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/misc/algorithm.hh"
#include "casm/monte/Conversions.hh"
#include "casm/monte/events/OccCandidate.hh"
#include "casm/monte/events/OccEventProposal.hh"
#include "casm/monte/events/OccLocation.hh"

namespace CASM {
namespace clexmonte {

namespace {
std::vector<monte::OccSwap>::const_iterator find_grand_canonical_swap(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_comp_n,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<Index> const &species_to_component_index_converter,
    MTRand &random_number_generator, monte::OccLocation const &occ_location,
    std::vector<monte::OccSwap>::const_iterator begin,
    std::vector<monte::OccSwap>::const_iterator end) {
  Eigen::VectorXd current_comp_n =
      composition_calculator.mean_num_each_component(occupation);
  auto const &index_converter = species_to_component_index_converter;

  double best_dist = (current_comp_n - target_comp_n).norm();

  double volume = occupation.size() / composition_calculator.n_sublat();
  double dn = 1. / volume;
  double tol = CASM::TOL;

  // store <distance_to_target_comp_n>:{swap_iterator, number of swaps}
  typedef std::vector<monte::OccSwap>::const_iterator iterator_type;
  typedef std::pair<iterator_type, Index> cand_and_count_pair;
  std::vector<cand_and_count_pair> choices;

  // check each possible swap for how close the composition is afterwards
  for (auto it = begin; it != end; ++it) {
    if (occ_location.cand_size(it->cand_a)) {
      Eigen::VectorXd tcomp_n = current_comp_n;
      tcomp_n[index_converter[it->cand_a.species_index]] -= dn;
      tcomp_n[index_converter[it->cand_b.species_index]] += dn;
      double dist = (tcomp_n - target_comp_n).norm();
      if (dist < best_dist - tol) {
        choices.clear();
        choices.push_back({it, occ_location.cand_size(it->cand_a)});
        best_dist = dist;
      } else if (dist < best_dist + tol) {
        choices.push_back({it, occ_location.cand_size(it->cand_a)});
      }
    }
  }
  if (!choices.size()) {
    return end;
  }

  // break ties randomly, weighted by number of candidates
  double sum = 0.0;
  for (const auto &val : choices) {
    sum += val.second;
  }

  double rand = random_number_generator.randExc(sum);
  sum = 0.0;
  for (const auto &val : choices) {
    sum += val.second;
    if (rand < sum) {
      return val.first;
    }
  }
  throw std::runtime_error(
      "Error in CASM::clexmonte::find_grand_canonical_swap, failed enforcing "
      "composition");
};

std::vector<Index> make_species_to_component_index_converter(
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

}  // namespace

/// \brief Apply grand canonical swaps to enforce target composition
///
/// Method:
/// - Find which of the provided grand canonical swap types transforms
///   the composition most closely to the target composition
/// - If no swap can improve the composition, return
/// - Propose and apply an event consistent with the found swap type
/// - Repeat
void enforce_composition(
    Eigen::VectorXi &occupation, Eigen::VectorXd const &target_comp_n,
    composition::CompositionCalculator const &composition_calculator,
    std::vector<monte::OccSwap> const &grand_canonical_swaps,
    monte::OccLocation &occ_location, MTRand &random_number_generator) {
  monte::Conversions const &convert = occ_location.convert();
  occ_location.initialize(occupation);

  // no guarantee convert species_index corresponds to comp_n index
  std::vector<Index> species_to_component_index_converter =
      make_species_to_component_index_converter(composition_calculator,
                                                convert);

  auto begin = grand_canonical_swaps.begin();
  auto end = grand_canonical_swaps.end();
  monte::OccEvent event;
  while (true) {
    auto it = find_grand_canonical_swap(
        occupation, target_comp_n, composition_calculator,
        species_to_component_index_converter, random_number_generator,
        occ_location, begin, end);

    if (it == end) {
      break;
    }

    /// propose event of chosen candidate type and apply swap
    monte::propose_grand_canonical_event_from_swap(event, occ_location, *it,
                                                   random_number_generator);
    occ_location.apply(event, occupation);
  }
}

}  // namespace clexmonte
}  // namespace CASM
