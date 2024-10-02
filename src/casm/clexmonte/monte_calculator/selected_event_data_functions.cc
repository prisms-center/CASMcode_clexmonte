#include "casm/clexmonte/monte_calculator/selected_event_data_functions.hh"

#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"

// debugging
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace clexmonte {
namespace monte_calculator {

namespace {
Index get_prim_event_index(
    std::shared_ptr<SelectedEvent> const &selected_event) {
  return selected_event->prim_event_data->prim_event_index;
}
std::vector<PrimEventData> const &get_prim_event_list(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return calculation->event_data().prim_event_list();
}
}  // namespace

/// \brief Make selected event count collecting function
/// ("selected_event.by_type")
monte::DiscreteVectorIntHistogramFunction make_selected_event_by_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // --- Data that will be used by the lamba function ---
  std::shared_ptr<std::vector<Index>> prim_event_index_to_index =
      std::make_shared<std::vector<Index>>();
  std::shared_ptr<SelectedEvent> selected_event = calculation->selected_event();
  std::shared_ptr<Eigen::VectorXi> event_type =
      std::make_shared<Eigen::VectorXi>(to_VectorXi(0));

  // Create:
  // - `key_to_index`: a map from event type name to event type index
  // - `value_labels`: a map from event type index (as a one element vector) to
  //    event type name (used for labeling the output)
  // - `prim_event_index_to_index`: a lookup table to convert from prim event
  //    index to event type index

  // get names in alphabetical order
  std::map<std::string, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[x.event_type_name] = 0;
  }

  // set index values
  std::map<Eigen::VectorXi, std::string, monte::LexicographicalCompare>
      value_labels;
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    value_labels.emplace(to_VectorXi(i_label), pair.first);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(key_to_index[x.event_type_name]);
  }

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_type",
      /* description */
      "Selected event count by event type",
      /* shape */
      {},
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event, event_type]() {
        Eigen::VectorXi &value = *event_type;
        value(0) =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return value;
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.by_equivalent_index")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_equivalent_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // --- Data that will be used by the lamba function ---
  std::shared_ptr<std::vector<Index>> prim_event_index_to_index =
      std::make_shared<std::vector<Index>>();
  std::shared_ptr<SelectedEvent> selected_event = calculation->selected_event();
  std::shared_ptr<Eigen::VectorXi> event_type =
      std::make_shared<Eigen::VectorXi>(to_VectorXi(0));

  // Create:
  // - `key_to_index`: a map from (event type name, equivalent_index) to event
  // type index
  // - `value_labels`: a map from event type index (as a one element vector) to
  //    event type name (used for labeling the output)
  // - `prim_event_index_to_index`: a lookup table to convert from prim event
  //    index to event type index

  // get names+equivalent_index in order
  std::map<std::pair<std::string, Index>, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[std::make_pair(x.event_type_name, x.equivalent_index)] = 0;
  }

  // set index values
  std::map<Eigen::VectorXi, std::string, monte::LexicographicalCompare>
      value_labels;
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    std::string label =
        pair.first.first + "." + std::to_string(pair.first.second);
    value_labels.emplace(to_VectorXi(i_label), label);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(
        key_to_index[std::make_pair(x.event_type_name, x.equivalent_index)]);
  }

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_equivalent_index",
      /* description */
      "Selected event count by event type and equivalent index",
      /* shape */
      {},
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event, event_type]() {
        Eigen::VectorXi &value = *event_type;
        value(0) =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return value;
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.by_prim_event_index")
monte::DiscreteVectorIntHistogramFunction
make_selected_event_by_prim_event_index_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // --- Data that will be used by the lamba function ---
  std::shared_ptr<std::vector<Index>> prim_event_index_to_index =
      std::make_shared<std::vector<Index>>();
  std::shared_ptr<SelectedEvent> selected_event = calculation->selected_event();
  std::shared_ptr<Eigen::VectorXi> event_type =
      std::make_shared<Eigen::VectorXi>(to_VectorXi(0));

  // Create:
  // - `key_to_index`: a map from (event type name, equivalent_index,
  //   is_forward) to event type index
  // - `value_labels`: a map from event type index (as a one element vector) to
  //    event type name (used for labeling the output)
  // - `prim_event_index_to_index`: a lookup table to convert from prim event
  //    index to event type index

  // get names+equivalent_index+is_forward in order
  std::map<std::tuple<std::string, Index, bool>, Index> key_to_index;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    key_to_index[std::make_tuple(x.event_type_name, x.equivalent_index,
                                 x.is_forward)] = 0;
  }

  // set index values
  std::map<Eigen::VectorXi, std::string, monte::LexicographicalCompare>
      value_labels;
  Index i_label = 0;
  for (auto &pair : key_to_index) {
    pair.second = i_label;
    std::string label = std::get<0>(pair.first) + "." +
                        std::to_string(std::get<1>(pair.first)) + "." +
                        (std::get<2>(pair.first) ? "forward" : "reverse");
    value_labels.emplace(to_VectorXi(i_label), label);
    ++i_label;
  }

  // create lookup table
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    prim_event_index_to_index->push_back(key_to_index[std::make_tuple(
        x.event_type_name, x.equivalent_index, x.is_forward)]);
  }

  // Create the function
  monte::DiscreteVectorIntHistogramFunction f(
      /* name */
      "selected_event.by_prim_event_index",
      /* description */
      "Selected event count by event type, equivalent index, and direction",
      /* shape */
      {},
      /* requires_event_state */
      false,
      /* function */
      [prim_event_index_to_index, selected_event, event_type]() {
        Eigen::VectorXi &value = *event_type;
        value(0) =
            (*prim_event_index_to_index)[get_prim_event_index(selected_event)];
        return value;
      },
      /* has_value_function */
      []() -> bool { return true; },
      /* max_size */
      10000);
  f.value_labels = value_labels;

  return f;
}

/// \brief Make selected event count collecting function
/// ("selected_event.<event_type>.by_equivalent_index")
std::vector<monte::DiscreteVectorIntHistogramFunction>
make_selected_event_by_equivalent_index_per_event_type_f(
    std::shared_ptr<MonteCalculator> const &calculation) {
  // prim_event_list should be present after calculation->reset():
  auto const &prim_event_list = get_prim_event_list(calculation);

  // --- Data that will be used by the lamba function ---
  std::shared_ptr<SelectedEvent> selected_event = calculation->selected_event();

  // Create:
  // - `key_to_index`: a map from event type name to event type index
  // - `value_labels`: a map from event type index (as a one element vector) to
  //    event type name (used for labeling the output)
  // - `prim_event_index_to_index`: a lookup table to convert from prim event
  //    index to event type index

  // get names in alphabetical order
  std::set<std::string> keys;
  for (clexmonte::PrimEventData const &x : prim_event_list) {
    keys.insert(x.event_type_name);
  }

  std::vector<monte::DiscreteVectorIntHistogramFunction> f_list;

  for (auto const &event_type_name : keys) {
    // --- Data that will be used by the lamba function ---
    auto value_ptr = std::make_shared<Eigen::VectorXi>(to_VectorXi(0));
    auto prim_event_index_to_index = std::make_shared<std::vector<Index>>();
    auto prim_event_index_to_has_value = std::make_shared<std::vector<bool>>();

    // create lookup tables and value labels
    std::map<Eigen::VectorXi, std::string, monte::LexicographicalCompare>
        value_labels;
    for (clexmonte::PrimEventData const &x : prim_event_list) {
      if (x.event_type_name == event_type_name) {
        prim_event_index_to_has_value->push_back(true);
        prim_event_index_to_index->push_back(x.equivalent_index);
        value_labels.emplace(
            to_VectorXi(x.equivalent_index),
            x.event_type_name + "." + std::to_string(x.equivalent_index));
      } else {
        prim_event_index_to_has_value->push_back(false);
        prim_event_index_to_index->push_back(-1);
      }
    }

    // Create the function
    monte::DiscreteVectorIntHistogramFunction f(
        /* name */
        "selected_event." + event_type_name + ".by_equivalent_index",
        /* description */
        "Selected event count by equivalent index for event_type " +
            event_type_name,
        /* shape */
        {},
        /* requires_event_state */
        false,
        /* function */
        [prim_event_index_to_index, selected_event, value_ptr]() {
          Eigen::VectorXi &value = *value_ptr;
          value(0) = (*prim_event_index_to_index)[get_prim_event_index(
              selected_event)];
          return value;
        },
        /* has_value_function */
        [prim_event_index_to_has_value, selected_event]() -> bool {
          return (*prim_event_index_to_has_value)[get_prim_event_index(
              selected_event)];
        },
        /* max_size */
        10000);
    f.value_labels = value_labels;
    f_list.push_back(f);
  }

  return f_list;
}

}  // namespace monte_calculator
}  // namespace clexmonte
}  // namespace CASM
