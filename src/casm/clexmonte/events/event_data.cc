#include "casm/clexmonte/events/event_data.hh"

#include "casm/clexmonte/misc/eigen.hh"

namespace CASM {
namespace clexmonte {

std::string name(PrimEventData const &prim_event_data) {
  std::string name = prim_event_data.event_type_name + "." +
                     std::to_string(prim_event_data.equivalent_index);
  if (prim_event_data.is_forward) {
    name += " (forward)";
  } else {
    name += " (reverse)";
  }
  return name;
}

SelectedEventInfo::SelectedEventInfo(
    std::vector<PrimEventData> const &_prim_event_list)
    : prim_event_list(_prim_event_list),
      prim_event_index_to_index(std::make_shared<std::vector<Index>>()),
      prim_event_index_to_has_value(std::make_shared<std::vector<bool>>()) {}

/// \brief Construct `prim_event_index_to_index` so that indices differentiate
///     by event type
void SelectedEventInfo::make_indices_by_type() {
  prim_event_index_to_index->clear();
  value_labels.clear();
  partition_names.clear();

  // get names in alphabetical order
  std::map<std::string, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[x.event_type_name] = 0;
  }

  // set index values for event type
  partition_names.resize(key_to_index.size());
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    partition_names[i_label] = pair.first;
    value_labels.emplace(to_VectorXl(i_label), pair.first);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(key_to_index[x.event_type_name]);
  }
}

/// \brief Construct `prim_event_index_to_index` so that indices differentiate
///     by event type and equivalent index
void SelectedEventInfo::make_indices_by_equivalent_index() {
  prim_event_index_to_index->clear();
  value_labels.clear();
  partition_names.clear();

  // get names+equivalent_index in order
  std::map<std::pair<std::string, Index>, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[std::make_pair(x.event_type_name, x.equivalent_index)] = 0;
  }

  // set index values
  partition_names.resize(key_to_index.size());
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    std::string label =
        pair.first.first + "." + std::to_string(pair.first.second);
    partition_names[i_label] = label;
    value_labels.emplace(to_VectorXl(i_label), label);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(
        key_to_index[std::make_pair(x.event_type_name, x.equivalent_index)]);
  }
}

/// \brief Construct `prim_event_index_to_index` so that indices differentiate
///     by event type, equivalent index, and direction
void SelectedEventInfo::make_indices_by_equivalent_index_and_direction() {
  prim_event_index_to_index->clear();
  value_labels.clear();
  partition_names.clear();

  // get names+equivalent_index+is_forward in order
  std::map<std::tuple<std::string, Index, bool>, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[std::make_tuple(x.event_type_name, x.equivalent_index,
                                 x.is_forward)] = 0;
  }

  // set index values
  partition_names.resize(key_to_index.size());
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    std::string label = std::get<0>(pair.first) + "." +
                        std::to_string(std::get<1>(pair.first)) + "." +
                        (std::get<2>(pair.first) ? "forward" : "reverse");
    partition_names[i_label] = label;
    value_labels.emplace(to_VectorXl(i_label), label);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(key_to_index[std::make_tuple(
        x.event_type_name, x.equivalent_index, x.is_forward)]);
  }
}

/// \brief Construct `prim_event_index_to_index` so that indices differentiate
///     by event type and equivalent index, but only for a single event type
///
/// - Does not make partition_names
/// - Other event types will have -1 for their index
void SelectedEventInfo::make_indices_by_equivalent_index_per_event_type(
    std::string event_type_name) {
  prim_event_index_to_index->clear();
  prim_event_index_to_has_value->clear();
  value_labels.clear();
  partition_names.clear();

  // create lookup tables and value labels
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    if (x.event_type_name == event_type_name) {
      prim_event_index_to_has_value->push_back(true);
      prim_event_index_to_index->push_back(x.equivalent_index);
      value_labels.emplace(
          to_VectorXl(x.equivalent_index),
          x.event_type_name + "." + std::to_string(x.equivalent_index));
    } else {
      prim_event_index_to_has_value->push_back(false);
      prim_event_index_to_index->push_back(-1);
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM