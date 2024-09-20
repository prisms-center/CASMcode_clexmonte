#ifndef CASM_clexmonte_MonteEventData
#define CASM_clexmonte_MonteEventData

#include "casm/clexmonte/monte_calculator/BaseMonteEventData.hh"

namespace CASM {
namespace clexmonte {

class MonteEventListIterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = EventID;
  using pointer = EventID const *;    // or also value_type*
  using reference = EventID const &;  // or also value_type&

  MonteEventListIterator() : m_data(nullptr), m_lib(nullptr), m_index(-1) {}

  MonteEventListIterator(std::shared_ptr<BaseMonteEventData> _data,
                         std::shared_ptr<RuntimeLibrary> _lib, bool _is_end)
      : m_data(_data), m_lib(_lib), m_index(m_data->new_iterator(_is_end)) {}

  ~MonteEventListIterator() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  MonteEventListIterator &operator=(MonteEventListIterator const &other) {
    m_data = other.m_data;
    m_lib = other.m_lib;
    m_index = m_data->copy_iterator(other.m_index);
    return *this;
  }

  reference operator*() const { return m_data->event_id(m_index); }
  pointer operator->() { return &m_data->event_id(m_index); }

  // Prefix increment (++it)
  MonteEventListIterator &operator++() {
    m_data->advance_iterator(m_index);
    return *this;
  }

  // Postfix increment (it++)
  MonteEventListIterator operator++(int) {
    MonteEventListIterator tmp = *this;
    tmp.m_index = m_data->copy_iterator(m_index);
    m_data->advance_iterator(m_index);
    return tmp;
  }

  bool operator==(MonteEventListIterator const &other) const {
    if (m_data == nullptr) {
      return false;
    }
    return m_data == other.m_data &&
           m_data->equal_iterator(m_index, other.m_index);
  }

  bool operator!=(MonteEventListIterator const &other) const {
    return !(*this == other);
  }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
  Index m_index;
};

/// Allows iterating over EventID
class MonteEventList {
 public:
  MonteEventList(std::shared_ptr<BaseMonteEventData> _data,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_data(_data), m_lib(_lib) {}

  ~MonteEventList() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  MonteEventListIterator begin() const {
    return MonteEventListIterator(m_data, m_lib, false);
  }

  MonteEventListIterator end() const {
    return MonteEventListIterator(m_data, m_lib, true);
  }

  /// The number of events
  Index size() const { return m_data->n_events(); }

  /// The current total event rate
  double total_rate() const { return m_data->total_rate(); }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
};

class MonteEventData {
 public:
  MonteEventData(std::shared_ptr<BaseMonteEventData> _data,
                 std::shared_ptr<RuntimeLibrary> _lib)
      : m_data(_data), m_lib(_lib), m_event_list(_data, _lib) {}

  ~MonteEventData() {
    // ensure BaseMonteEventData is deleted before library
    m_data.reset();
  }

  /// The system
  std::shared_ptr<system_type> system() const { return m_data->system; }

  /// The `prim events`, one translationally distinct instance
  /// of each event, associated with origin primitive cell
  std::vector<clexmonte::PrimEventData> const &prim_event_list() const {
    return m_data->prim_event_list;
  }

  /// Information about what sites may impact each prim event
  std::vector<clexmonte::EventImpactInfo> const &prim_impact_info_list() const {
    return m_data->prim_impact_info_list;
  }

  /// Get the formation energy coefficients
  clexulator::SparseCoefficients const &formation_energy_coefficients() const {
    return m_data->formation_energy_coefficients();
  }

  /// Get the attempt frequency coefficients for a specific event
  clexulator::SparseCoefficients const &freq_coefficients(
      Index prim_event_index) const {
    return m_data->freq_coefficients(prim_event_index);
  }

  /// Get the KRA coefficients for a specific event
  clexulator::SparseCoefficients const &kra_coefficients(
      Index prim_event_index) const {
    return m_data->kra_coefficients(prim_event_index);
  }

  // -- Event list --

  /// Get EventID
  MonteEventList const &event_list() const { return m_event_list; }

  // -- Event info (accessed by EventID) --

  /// The monte::OccEvent that can apply the event for the current state of the
  /// internal iterator
  monte::OccEvent const &event_to_apply(EventID const &id) const {
    return m_data->event_to_apply(id);
  }

  /// Return the current rate for a specific event
  double event_rate(EventID const &id) const { return m_data->event_rate(id); }

  /// Calculate event state data
  EventState const &event_state(EventID const &id) const {
    return m_data->event_state(id);
  }

  /// The events that must be updated if the specified event occurs
  std::vector<EventID> const &event_impact(EventID const &id) const {
    return m_data->impact(id);
  }

 private:
  std::shared_ptr<BaseMonteEventData> m_data;
  std::shared_ptr<RuntimeLibrary> m_lib;
  MonteEventList m_event_list;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
