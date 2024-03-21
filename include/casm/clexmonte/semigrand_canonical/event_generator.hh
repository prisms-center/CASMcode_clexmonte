#ifndef CASM_clexmonte_semigrand_canonical_event_generator
#define CASM_clexmonte_semigrand_canonical_event_generator

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/semigrand_canonical/conditions.hh"
#include "casm/clexmonte/semigrand_canonical/potential.hh"
#include "casm/monte/RandomNumberGenerator.hh"
#include "casm/monte/events/OccEventProposal.hh"

namespace CASM {
namespace clexmonte {
namespace semigrand_canonical {

/// \brief Propose and apply semi-grand canonical events
template <typename EngineType>
class SemiGrandCanonicalEventGenerator {
 public:
  /// \brief Constructor
  ///
  /// Notes:
  /// - One and only one of `_semigrand_canonical_swaps` and
  /// `_semigrand_canonical_multiswaps` should have size != 0
  ///
  /// \param _semigrand_canonical_swaps Single site swap types for semi-grand
  ///     canonical Monte Carlo events. If size > 0, only these events are
  ///     proposed.
  /// \param _semigrand_canonical_multiswaps Multiple site swap types for
  ///     semi-grand canonical Monte Carlo events, such as charge neutral
  ///     events. These events are only proposed if no single swaps are
  ///     provided.
  SemiGrandCanonicalEventGenerator(
      std::vector<monte::OccSwap> const &_semigrand_canonical_swaps,
      std::vector<monte::MultiOccSwap> const &_semigrand_canonical_multiswaps)
      : state(nullptr),
        occ_location(nullptr),
        semigrand_canonical_swaps(_semigrand_canonical_swaps),
        semigrand_canonical_multiswaps(_semigrand_canonical_multiswaps),
        use_multiswaps(semigrand_canonical_swaps.size() == 0) {
    if (semigrand_canonical_swaps.size() == 0 &&
        semigrand_canonical_multiswaps.size() == 0) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalEventGenerator: "
          "semigrand_canonical_swaps.size() == 0 && "
          "semigrand_canonical_multiswaps.size() == 0");
    }
    if (semigrand_canonical_swaps.size() != 0 &&
        semigrand_canonical_multiswaps.size() != 0) {
      throw std::runtime_error(
          "Error in SemiGrandCanonicalEventGenerator: "
          "semigrand_canonical_swaps.size() != 0 && "
          "semigrand_canonical_multiswaps.size() != 0");
    }
  }

  /// \brief The current state for which events are proposed and applied. Can be
  ///     nullptr, but must be set for use.
  state_type *state;

  /// Occupant tracker
  monte::OccLocation *occ_location;

  /// \brief Single swap types for semi-grand canonical Monte Carlo events
  std::vector<monte::OccSwap> semigrand_canonical_swaps;

  /// \brief Multiple swap types for semi-grand canonical Monte Carlo events
  std::vector<monte::MultiOccSwap> semigrand_canonical_multiswaps;

  /// \brief If true, propose events from multiswaps, else propose events from
  /// single swaps
  bool use_multiswaps;

  /// \brief The current proposed event
  monte::OccEvent occ_event;

 public:
  /// \brief Set the current Monte Carlo state and occupant locations
  ///
  /// Notes:
  /// - Must be called before `propose` or `apply`
  ///
  /// \param _state The current state for which events are proposed and applied.
  ///     Throws if nullptr.
  /// \param _occ_location An occupant location tracker, which enables efficient
  ///     event proposal. It must already be initialized with the input state.
  ///     Throws if nullptr.
  void set(state_type *_state, monte::OccLocation *_occ_location) {
    this->state =
        throw_if_null(_state,
                      "Error in SemiGrandCanonicalEventGenerator::set: "
                      "_state==nullptr");
    this->occ_location =
        throw_if_null(_occ_location,
                      "Error in SemiGrandCanonicalEventGenerator::set: "
                      "_occ_location==nullptr");
  }

  /// \brief Propose a Monte Carlo occupation event, returning a reference
  ///
  /// Notes:
  /// - Must call `set` before `propose` or `apply`
  ///
  /// \param random_number_generator A random number generator
  monte::OccEvent const &propose(
      monte::RandomNumberGenerator<EngineType> &random_number_generator) {
    if (this->use_multiswaps) {
      return monte::propose_semigrand_canonical_multiswap_event(
          this->occ_event, *this->occ_location,
          this->semigrand_canonical_multiswaps, random_number_generator);
    } else {
      return monte::propose_semigrand_canonical_event(
          this->occ_event, *this->occ_location, this->semigrand_canonical_swaps,
          random_number_generator);
    }
  }

  /// \brief Update the occupation of the current state using the provided event
  void apply(monte::OccEvent const &e) {
    this->occ_location->apply(e, get_occupation(*this->state));
  }
};

}  // namespace semigrand_canonical
}  // namespace clexmonte
}  // namespace CASM

#endif
