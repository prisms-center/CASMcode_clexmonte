#ifndef CASM_clexmonte_state_parse_conditions_impl
#define CASM_clexmonte_state_parse_conditions_impl

#include <memory>

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/state/io/json/parse_conditions.hh"
#include "casm/clexmonte/state/make_conditions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/monte/ValueMap.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace clexmonte {

/// \brief Parse 'temperature' (as a required attribute)
///
/// If successfully parsed, `parser->value` will contain a
/// monte::ValueMap with:
/// - scalar_values["temperature"]: (size 1)
///
/// If unsuccesfully parsed, `parser.valid() == false`.
///
/// Expected:
///
///   "temperature": number (optional)
///     Temperature in K.
///
/// \tparam ConditionsType The ConditionsType must support
///     `ConditionsType::set_temperature(double temperature)`.
/// \param parser InputParser, which must have non-empty value
///
template <typename ConditionsType>
void parse_temperature(InputParser<ConditionsType> &parser) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_temperature: parser must have non-empty value");
  }
  auto subparser = parser.template subparse<double>("temperature");
  if (subparser->valid()) {
    parser.value->set_temperature(*subparser->value);
  }
}

/// \brief Parse 'param_chem_pot' (as a required attribute)
///
/// \tparam ConditionsType The ConditionsType must support
///     `ConditionsType::set_param_composition(Eigen::VectorXd
///     param_composition, bool is_increment)`.
/// \param parser Conditions parser
/// \param tol Tolerance used to check equivalence if 'mol_composition' and
///     'param_composition' are both included
///
/// Expected:
///
///   "param_chem_pot": Union[list,dict,None]
///     Potential conjugate to parametric composition, in terms of the
///     chosen composition axes. The keys are the axes names ("a", "b",
///     etc.), and values are the corresponding potential value. If an
///     array, must match composition axes size. All composition axes
///     must be included.
///
template <typename ConditionsType>
void parse_param_chem_pot(
    InputParser<ConditionsType> &parser,
    composition::CompositionConverter const &composition_converter) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_composition: parser must have non-empty value");
  }

  if (!parser.self.contains("param_chem_pot")) {
    std::stringstream msg;
    msg << "Error: 'param_chem_pot' must be provided";
    parser.insert_error("param_chem_pot", msg.str());
    return;
  }

  Eigen::VectorXd param_chem_pot;

  if (parser.self["param_chem_pot"].is_array()) {
    parser.optional(param_chem_pot, "param_chem_pot");
    if (param_chem_pot.size() !=
        composition_converter.independent_compositions()) {
      std::stringstream msg;
      msg << "Error: 'param_chem_pot' size mismatch.";
      parser.insert_error("param_chem_pot", msg.str());
    }
  } else if (parser.self["param_chem_pot"].is_object()) {
    std::map<std::string, double> input;
    parser.optional(input, "param_chem_pot");
    try {
      param_chem_pot = make_param_chem_pot(composition_converter, input);
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error: could not construct composition from option "
             "'param_chem_pot'.";
      parser.insert_error("param_chem_pot", e.what());
    }
  } else {
    std::stringstream msg;
    msg << "Error: 'param_chem_pot' must be an array or object";
    parser.insert_error("param_chem_pot", msg.str());
  }

  if (parser.valid()) {
    parser.value->set_param_chem_pot(param_chem_pot);
  }
}

/// \brief Parse 'mol_composition' and/or 'param_composition' (at least one of
/// which is required, if both are provided they must be equivalent)
///
/// \tparam ConditionsType The ConditionsType must support
///     `ConditionsType::set_param_composition(Eigen::VectorXd
///     param_composition, bool is_increment)`.
/// \param parser
/// \param is_increment
/// \param tol Tolerance used to check equivalence if 'mol_composition' and
///     'param_composition' are both included
///
/// Expected:
///
///   "mol_composition": array or dict (optional)
///     Composition in number per primitive cell. If a dict, the keys are the
///     component names, and values are the number of that component per
///     primitive cell. If an array, must match composition_converter components
///     order and size. All components in the system must be included. Must sum
///     to the number of sites per prim cell.
///
///   "param_composition": Union[list,dict,None]
///     Parametric composition, in terms of the chosen composition axes. The
///     keys are the axes names ("a", "b", etc.), and values are the
///     corresponding parametric composition value. If an array, must match
///     composition axes size. All composition axes must be included.
///
template <typename ConditionsType>
void parse_composition(
    InputParser<ConditionsType> &parser,
    composition::CompositionConverter const &composition_converter,
    bool is_increment) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_composition: parser must have non-empty value");
  }

  Eigen::VectorXd tmp_vectorxd;
  std::optional<Eigen::VectorXd> mol_composition;
  std::optional<Eigen::VectorXd> param_composition;

  /// Get "mol_composition" if exists
  if (parser.self.contains("mol_composition")) {
    if (parser.self["mol_composition"].is_array()) {
      parser.optional(tmp_vectorxd, "mol_composition");
      mol_composition = tmp_vectorxd;
    } else if (parser.self["mol_composition"].is_obj()) {
      try {
        std::map<std::string, double> input;
        parser.optional(input, "mol_composition");
        bool do_not_convert = true;
        if (!is_increment) {
          tmp_vectorxd = make_mol_composition(composition_converter, input,
                                              do_not_convert);
        } else {
          tmp_vectorxd = make_mol_composition_increment(composition_converter,
                                                        input, do_not_convert);
        }
        mol_composition = tmp_vectorxd;
      } catch (std::exception &e) {
        std::stringstream msg;
        msg << "Error: could not construct composition from option "
               "'mol_composition'. "
            << e.what();
        parser.insert_error("mol_composition", msg.str());
        return;
      }
    } else {
      std::stringstream msg;
      msg << "Error: 'mol_composition' must be an array or object";
      parser.insert_error("mol_composition", msg.str());
      return;
    }
  }

  /// Get "param_composition" if exists
  if (parser.self.contains("param_composition")) {
    if (parser.self["param_composition"].is_array()) {
      parser.optional(tmp_vectorxd, "param_composition");
      param_composition = tmp_vectorxd;
    } else if (parser.self["param_composition"].is_object()) {
      std::map<std::string, double> input;
      parser.optional(input, "param_composition");
      try {
        bool do_not_convert = true;
        if (!is_increment) {
          tmp_vectorxd = make_param_composition(composition_converter, input,
                                                do_not_convert);
        } else {
          tmp_vectorxd = make_param_composition_increment(
              composition_converter, input, do_not_convert);
        }
        param_composition = tmp_vectorxd;
      } catch (std::exception &e) {
        std::stringstream msg;
        msg << "Error: could not construct composition from option "
               "'param_composition'.";
        parser.insert_error("param_composition", e.what());
        return;
      }
    } else {
      std::stringstream msg;
      msg << "Error: 'param_composition' must be an array or object";
      parser.insert_error("param_composition", msg.str());
      return;
    }
  }

  if (parser.valid()) {
    parser.value->set_composition(param_composition, mol_composition,
                                  is_increment);
  }
}

}  // namespace clexmonte
}  // namespace CASM

#endif
