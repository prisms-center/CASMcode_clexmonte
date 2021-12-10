#include "casm/clexmonte/system/io/json/parse_composition.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace clexmonte {

namespace {

/// \brief Parse comp_n from an array, validating size and types
///
/// Expects `parser.self` has form such as:
/// \code
/// {
///   "comp_n": [2.0, 1.6, 0.4],
///   ...
/// }
/// \endcode
///
/// Requires:
/// - size of "comp_n" matches composition_converter
void parse_comp_n_array(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter) {
  std::unique_ptr<Eigen::VectorXd> comp_n =
      parser.optional<Eigen::VectorXd>("comp_n");
  if (comp_n == nullptr) {
    return;
  }
  Index n_components = composition_converter.components().size();
  if (comp_n->size() != n_components) {
    parser.insert_error("comp_n",
                        "Error parsing \"comp_n\": size != the number "
                        "components in this system.");
    return;
  }

  parser.value = std::move(comp_n);
}

/// \brief Parse comp_x from an array, validating size and types
///
/// Expects `parser.self` has form such as:
/// \code
/// {
///   "comp_x": [0.5, 0.3],
///   ...
/// }
/// \endcode
///
/// Requires:
/// - size of "comp_x" matches composition_converter
void parse_comp_x_array(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter) {
  std::unique_ptr<Eigen::VectorXd> comp_x =
      parser.optional<Eigen::VectorXd>("comp_x");
  if (comp_x == nullptr) {
    return;
  }
  Index n_axes = composition_converter.independent_compositions();
  if (comp_x->size() != n_axes) {
    parser.insert_error("comp_x",
                        "Error parsing \"comp_x\": size != the number of "
                        "independent composition axes in this system.");
    return;
  }

  parser.value = std::move(comp_x);
}

/// \brief Parse comp_n from an object, validating size and types
///
/// Expects `parser.self` has form such as:
/// \code
/// {
///   "comp_n": {
///      "Zr": 2.0,
///      "O: 1.6,
///      "Va": 0.4
///   }
/// }
/// \endcode
///
/// Requires:
/// - size of "comp_n" matches composition_converter
void parse_comp_n_object(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter) {
  // parse comp_n object
  // example: "comp_n": {"A": 1.5, "B": 0.5}
  std::vector<std::string> components = composition_converter.components();

  typedef std::map<std::string, double> expected_type;
  std::unique_ptr<expected_type> comp_n_map =
      parser.optional<expected_type>("comp_n");
  if (comp_n_map == nullptr) {
    return;
  }
  Index n_components = components.size();
  if (comp_n_map->size() != n_components) {
    parser.insert_error("comp_n",
                        "Error parsing \"comp_n\": size != the number "
                        "components in this system.");
    return;
  }

  // check for missing components
  bool valid = true;
  for (auto const &name : components) {
    auto map_it = comp_n_map->find(name);
    auto map_end = comp_n_map->end();
    if (map_it == map_end) {
      std::stringstream msg;
      msg << "Error parsing \"comp_n\": missing component \"" << name << "\"";
      parser.insert_error("comp_n", msg.str());
      valid = false;
    }
  }

  // read component compositions
  auto vector_begin = components.begin();
  auto vector_end = components.end();
  Eigen::VectorXd comp_n(n_components);
  for (auto const &pair : *comp_n_map) {
    std::string name = pair.first;
    double value = pair.second;
    auto vector_it = std::find(vector_begin, vector_end, name);
    if (vector_it == vector_end) {
      std::stringstream msg;
      msg << "Error parsing \"comp_n\": unexpected component \"" << name
          << "\"";
      parser.insert_error("comp_n", msg.str());
      valid = false;
      continue;
    }
    comp_n(std::distance(vector_begin, vector_it)) = value;
  }

  if (!valid) {
    return;
  }

  parser.value = std::make_unique<Eigen::VectorXd>(comp_n);
}

/// \brief Parse comp_x from an object, validating size and types
///
/// Expects `parser.self` has form such as:
/// \code
/// {
///   "comp_x": {
///      "a": 0.5,
///      "b: 0.3
///   }
/// }
/// \endcode
///
/// Requires:
/// - size of "comp_x" matches composition_converter
void parse_comp_x_object(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter) {
  typedef std::map<std::string, double> expected_type;
  std::unique_ptr<expected_type> comp_x_map =
      parser.optional<expected_type>("comp_x");
  if (comp_x_map == nullptr) {
    return;
  }
  Index n_axes = composition_converter.independent_compositions();
  if (comp_x_map->size() != n_axes) {
    parser.insert_error("comp_x",
                        "Error parsing \"comp_x\": size != the number of "
                        "independent composition axes in this system.");
    return;
  }

  // check for missing axes
  bool valid = true;
  std::vector<std::string> axis_names;
  for (Index i = 0; i < n_axes; ++i) {
    std::string name = composition_converter.comp_var(i);
    axis_names.push_back(name);
    auto map_it = comp_x_map->find(name);
    auto map_end = comp_x_map->end();
    if (map_it == map_end) {
      std::stringstream msg;
      msg << "Error parsing \"comp_x\": missing composition for axis \"" << name
          << "\"";
      parser.insert_error("comp_x", msg.str());
      valid = false;
    }
  }

  // read component compositions
  auto vector_begin = axis_names.begin();
  auto vector_end = axis_names.end();
  Eigen::VectorXd comp_x(n_axes);
  for (auto const &pair : *comp_x_map) {
    std::string name = pair.first;
    double value = pair.second;
    auto vector_it = std::find(vector_begin, vector_end, name);
    if (vector_it == vector_end) {
      std::stringstream msg;
      msg << "Error parsing \"comp_x\": unexpected name \"" << name << "\"";
      parser.insert_error("comp_x", msg.str());
      valid = false;
      continue;
    }
    comp_x(std::distance(vector_begin, vector_it)) = value;
  }

  if (!valid) {
    return;
  }

  parser.value = std::make_unique<Eigen::VectorXd>(comp_x);
}

}  // namespace

/// \brief Parse composition (either comp_n or comp_x) from JSON
///
/// Expects input of form documented by `parse_composition_for_comp_n`
///
/// Parses composition from JSON, output as a vector. If successfully parsed,
/// `parser->value` will contain an Eigen::VectorXd.
///
/// \param parser JSON InputParser. If successfully parsed, `parser->value` will
///     contain an `Eigen::VectorXd`.
/// \param composition_converter Composition axes and conversions
/// \param is_comp_n Will be set to true if successfully parsed from "comp_n".
///     Will be set to false if successfully parsed from "comp_x".
///
void parse_composition_as_provided(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool &is_comp_n) {
  if (parser.self.contains("comp_n") == parser.self.contains("comp_x")) {
    parser.insert_error(
        "comp_n", "Error: requires exactly one of \"comp_n\", \"comp_x\"");
    return;
  }

  auto it = parser.self.find("comp_n");
  auto end = parser.self.end();
  if (it != end) {
    if (it->is_array()) {
      parse_comp_n_array(parser, composition_converter);
    } else if (it->is_object()) {
      parse_comp_n_object(parser, composition_converter);
    } else {
      parser.insert_error("comp_n",
                          "Error parsing \"comp_n\": must be array or object");
      return;
    }
    if (parser.value != nullptr) {
      is_comp_n = true;
    }
    return;
  }

  it = parser.self.find("comp_x");
  if (it != end) {
    if (it->is_array()) {
      parse_comp_x_array(parser, composition_converter);
    } else if (it->is_object()) {
      parse_comp_x_object(parser, composition_converter);
    } else {
      parser.insert_error("comp_n",
                          "Error parsing \"comp_n\": must be array or object");
      return;
    }
    if (parser.value != nullptr) {
      is_comp_n = false;
    }
    return;
  }
}

/// \brief Parse composition from JSON, return as comp_n or comp_n increment
///
/// \param parser JSON InputParser. If successfully parsed, `parser->value` will
///     contain an `Eigen::VectorXd` represent comp_x or comp_x increment.
/// \param composition_converter Composition axes and conversions
/// \param is_increment If true, parse input as comp_n increment. Otherwise,
///     parse as comp_n.
///
/// Expected JSON format:
///
///   "comp_n": array of number or dict (optional)
///     Composition in number per primitive cell. May be:
///
///     - An array of number, specifying the number of each component per
///       primitive cell, interpreted using the order of the `"components"`
///       specified for composition axes. The size must match the number of
///       components.
///
///     - A dict, where the keys are the component names, and values are
///       the number of that component per primitive cell. All components in
///       the system must be included.
///
///   "comp_x": array of number or dict (optional)
///     Parametric composition, in terms of the chosen composition axes. Will
///     be converted to `"comp_n"`. May
///     be:
///
///     - An array of number, corresponding to
///       `[comp_x(a), comp_x(b), ...]`. The size must match the number of
///       composition axes.
///
///     - A dict, where the keys are the axes names ("a", "b", etc.), and
///       values are the corresponding parametric composition value.
///       All composition axes must be included.
///
void parse_composition_for_comp_n(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment) {
  bool is_comp_n;
  parse_composition_as_provided(parser, composition_converter, is_comp_n);

  if (parser.value == nullptr) {
    return;
  }

  if (is_comp_n) {
    double expected_comp_n_sum =
        is_increment ? 0.0 : composition_converter.origin().sum();
    // check sum is expected value
    // (depending on context, should be number of sublatticese or zero)
    if (!CASM::almost_equal(parser.value->sum(), expected_comp_n_sum, TOL)) {
      std::stringstream msg;
      msg << "Error parsing \"comp_n\": sum != " << expected_comp_n_sum;
      parser.insert_error("comp_n", msg.str());
      parser.value.reset();
      return;
    }
  } else {
    if (!is_increment) {
      *parser.value = composition_converter.mol_composition(*parser.value);
    } else {
      *parser.value = composition_converter.dmol_composition(*parser.value);
    }
  }
}

/// \brief Parse composition from JSON, return as comp_x or comp_x increment
///
/// Expects input of form documented by `parse_composition_for_comp_n`
///
/// \param parser JSON InputParser. If successfully parsed, `parser->value` will
///     contain an `Eigen::VectorXd` represent comp_x or comp_x increment.
/// \param composition_converter Composition axes and conversions
/// \param is_increment If true, parse input as comp_x increment. Otherwise,
///     parse as comp_x.
///
void parse_composition_for_comp_x(
    InputParser<Eigen::VectorXd> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment) {
  bool is_comp_n;
  parse_composition_as_provided(parser, composition_converter, is_comp_n);

  if (parser.value == nullptr) {
    return;
  }

  if (is_comp_n) {
    double expected_comp_n_sum =
        is_increment ? 0.0 : composition_converter.origin().sum();
    // check sum is expected value
    // (depending on context, should be number of sublatticese or zero)
    if (!CASM::almost_equal(parser.value->sum(), expected_comp_n_sum, TOL)) {
      std::stringstream msg;
      msg << "Error parsing \"comp_n\": sum != " << expected_comp_n_sum;
      parser.insert_error("comp_n", msg.str());
      parser.value.reset();
      return;
    }
    if (!is_increment) {
      *parser.value = composition_converter.param_composition(*parser.value);
    } else {
      *parser.value = composition_converter.dparam_composition(*parser.value);
    }
  }
}

}  // namespace clexmonte
}  // namespace CASM
