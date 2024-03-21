#include "casm/clexmonte/state/io/json/parse_conditions.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/misc/eigen.hh"
#include "casm/clexmonte/state/Conditions.hh"
#include "casm/clexmonte/state/CorrMatchingPotential.hh"
#include "casm/clexmonte/state/io/json/CorrMatchingPotential_json_io.hh"
#include "casm/clexmonte/state/make_conditions.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/composition/CompositionConverter.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/monte/ValueMap.hh"

namespace CASM {
namespace clexmonte {

void parse(InputParser<Conditions> &parser,
           std::shared_ptr<system_type> const &system, bool is_increment) {
  auto conditions_subparser = parser.parse_as_with<monte::ValueMap>(
      parse_conditions, system, is_increment);
  if (conditions_subparser->valid()) {
    if (!is_increment) {
      parser.value =
          std::make_unique<Conditions>(make_conditions_from_value_map(
              *conditions_subparser->value, *get_prim_basicstructure(*system),
              get_composition_converter(*system),
              get_random_alloy_corr_f(*system), CASM::TOL /*TODO*/));
    } else {
      parser.value =
          std::make_unique<Conditions>(make_conditions_increment_from_value_map(
              *conditions_subparser->value, *get_prim_basicstructure(*system),
              get_composition_converter(*system),
              get_random_alloy_corr_f(*system), CASM::TOL /*TODO*/));
    }
  }
}

/// \brief Parse all conditions
///
/// Example input (all optional for purposes of this method):
////// \code
/// {
///   "temperature" : number,
///
///   // specify parametric composition, referenced to composition axes of
///   composition_converter
///   // option 1: object, specifying value with axis name for key
///   "param_composition" : {"a": number, "b": number, ...},
///   // option 1: array, specifying value for ["a", "b", ...]
///   "param_composition" : [number, number, ...],
///
///   // specify mol composition (number per primitive cell), using one of two
///   options:
///   // option 1: object, specifying values with key from
///   // composition_converter.components()
///   "mol_composition" : {"A": number, "B": number, ...},
///   // option 2: array, in order of composition_converter.components()
///   "mol_composition" : [number, number, ...],
///
///   // specify potential conjugate to composition axes
///   // option 1: object, specifying value with axis name for key
///   "param_chem_pot" : {"a": number, "b": number, ...},
///   // option 1: array, specifying value for ["a", "b", ...]
///   "param_chem_pot" : [number, number, ...],
///
///   // For use with correlation-matching potentials or other custom
///   potentials,
///   // the formation energy contribution to the potential can be turned off.
///   // The default value is true, except default=false if one of the following
///   // are included:
///   // - "corr_matching_pot"
///   // - "random_alloy_corr_matching_pot"
///   "include_formation_energy": bool,
///
///   // quadratic potential of parameteric composition, minimum location
///   "param_comp_quad_pot_target": [number, number, ...],
///
///   // potential curvature (vector, for on-diagonal terms only)
///   "param_comp_quad_pot_vector": [number, number, ...]
///
///   // potential curvature (matrix, to include off-diagonal terms)
///   "param_comp_quad_pot_matrix": [
///     [number, number, ...],
///     [number, number, ...],
///     ...],
///
///   // linear potential of order parameter
///   "order_parameter_pot": [number, number, ...],
///
///   // quadratic potential of order parameter, minimum location
///   "order_parameter_quad_pot_target": [number, number, ...],
///
///   // potential curvature (vector, for on-diagonal terms only)
///   "order_parameter_quad_pot_vector": [number, number, ...]
///
///   // potential curvature (matrix, to include off-diagonal terms)
///   "order_parameter_quad_pot_matrix": [ // curvature (matrix)
///     [number, number, ...],
///     [number, number, ...],
///     ...],
///
///   // correlation-matching potential (SQS generation)
///   "corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "targets": [
///       {"index":int, "value":number, "weight":number (default=1.)},
///       ...
///     ]
///   }
///
///   // random-alloy correlation-matching potential (SQS generation)
///   "random_alloy_corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "sublattice_prob": [
///       [number, number, ...], // sublattice 0
///       [number, number, ...], // sublattice 1
///       ...
///     ],
///     "targets": [
///       {"index":int, "weight":number (default=1.)},
///       ...
///     ]
///   }
/// }
/// \endcode
///
void parse_conditions(InputParser<monte::ValueMap> &parser,
                      std::shared_ptr<system_type> const &system,
                      bool is_increment) {
  if (parser.value == nullptr) {
    parser.value = std::make_unique<monte::ValueMap>();
  }
  parse_temperature(parser);
  parse_mol_composition(parser, system, is_increment);
  parse_param_composition(parser, system, is_increment);
  parse_param_chem_pot(parser, system);

  parse_include_formation_energy(parser);

  parse_vector(parser, "param_comp_quad_pot_target");
  parse_vector(parser, "param_comp_quad_pot_vector");
  parse_matrix(parser, "param_comp_quad_pot_matrix");

  parse_vector(parser, "order_parameter_pot");

  parse_vector(parser, "order_parameter_quad_pot_target");
  parse_vector(parser, "order_parameter_quad_pot_vector");
  parse_matrix(parser, "order_parameter_quad_pot_matrix");

  parse_corr_matching_pot(parser, is_increment);
  parse_random_alloy_corr_matching_pot(parser, system, is_increment);

  if (!parser.valid()) {
    return;
  }

  try {
    // This is an awkward way to do things... but it allows input of
    // mol_composition or param_composition and checks consistency if both exist
    if (!is_increment) {
      Conditions conditions = make_conditions_from_value_map(
          *parser.value, *get_prim_basicstructure(*system),
          get_composition_converter(*system), get_random_alloy_corr_f(*system),
          CASM::TOL /*TODO*/);
      *parser.value = make_value_map_from_conditions(conditions);
    } else {
      Conditions conditions_increment =
          make_conditions_increment_from_value_map(
              *parser.value, *get_prim_basicstructure(*system),
              get_composition_converter(*system),
              get_random_alloy_corr_f(*system), CASM::TOL /*TODO*/);
      *parser.value =
          make_value_map_from_conditions_increment(conditions_increment);
    }
  } catch (std::exception &e) {
    parser.error.insert(e.what());
    parser.value.reset();
  }
}

/// \brief Parse temperature scalar value
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
/// \param parser InputParser, which must have non-empty value
///
void parse_temperature(InputParser<monte::ValueMap> &parser) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_temperature: parser must have non-empty value");
  }
  parser.optional(parser.value->scalar_values["temperature"], "temperature");
}

/// \brief Parse "mol_composition" and store as
///     "mol_composition" vector values
///
/// If successfully parsed, `parser->value` will contain a
/// monte::ValueMap with:
/// - vector_values["mol_composition"]: (size = system components size)
///
/// If unsuccessfully parsed, `parser.valid() == false`.
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
/// Requires:
/// - get_composition_converter(system_type const &system);
void parse_mol_composition(InputParser<monte::ValueMap> &parser,
                           std::shared_ptr<system_type> const &system,
                           bool is_increment) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_mol_composition: parser must have non-empty value");
  }

  auto const &composition_converter = get_composition_converter(*system);
  if (parser.self["mol_composition"].is_array()) {
    Eigen::VectorXd &mol_composition =
        parser.value->vector_values["mol_composition"];
    parser.optional(mol_composition, "mol_composition");
    if (mol_composition.size() != composition_converter.components().size()) {
      std::stringstream msg;
      msg << "Error: 'mol_composition' size mismatch.";
      parser.insert_error("mol_composition", msg.str());
    }
  } else if (parser.self["mol_composition"].is_obj()) {
    std::map<std::string, double> input;
    parser.optional(input, "mol_composition");
    try {
      bool do_not_convert = true;
      if (!is_increment) {
        parser.value->vector_values["mol_composition"] =
            make_mol_composition(composition_converter, input, do_not_convert);
      } else {
        parser.value->vector_values["mol_composition"] =
            make_mol_composition_increment(composition_converter, input,
                                           do_not_convert);
      }
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error: could not construct composition from option "
             "'mol_composition'. "
          << e.what();
      parser.insert_error("mol_composition", msg.str());
    }
  } else {
    std::stringstream msg;
    msg << "Error: 'mol_composition' must be an array or object";
    parser.insert_error("mol_composition", msg.str());
  }
}

/// \brief Parse a vector or dict of parametric composition : value
///
/// If successfully parsed, `parser->value` will contain a
/// monte::ValueMap with:
/// - vector_values["param_composition"]: (size = system components size)
///
/// If unsuccessfully parsed, `parser.valid() == false`.
///
/// Expected:
///
///   "param_composition": Union[list,dict,None]
///     Parametric composition, in terms of the chosen composition axes. The
///     keys are the axes names ("a", "b", etc.), and values are the
///     corresponding parametric composition value. If an array, must match
///     composition axes size. All composition axes must be included.
///
/// Requires:
/// - get_composition_converter(system_type const &system);
void parse_param_composition(InputParser<monte::ValueMap> &parser,
                             std::shared_ptr<system_type> const &system,
                             bool is_increment) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_param_composition: parser must have non-empty value");
  }

  if (!parser.self.contains("param_composition")) {
    return;
  }

  auto const &composition_converter = get_composition_converter(*system);
  if (parser.self["param_composition"].is_array()) {
    Eigen::VectorXd &value = parser.value->vector_values["param_composition"];
    parser.optional(value, "param_composition");
    if (value.size() != composition_converter.independent_compositions()) {
      std::stringstream msg;
      msg << "Error: 'param_composition' size mismatch.";
      parser.insert_error("param_composition", msg.str());
    }
  } else if (parser.self["param_composition"].is_object()) {
    std::map<std::string, double> input;
    parser.optional(input, "param_composition");
    try {
      bool do_not_convert = true;
      if (!is_increment) {
        parser.value->vector_values["param_composition"] =
            make_param_composition(composition_converter, input,
                                   do_not_convert);
      } else {
        parser.value->vector_values["param_composition"] =
            make_param_composition_increment(composition_converter, input,
                                             do_not_convert);
      }
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error: could not construct composition from option "
             "'param_composition'.";
      parser.insert_error("param_composition", e.what());
    }
  } else {
    std::stringstream msg;
    msg << "Error: 'param_composition' must be an array or object";
    parser.insert_error("param_composition", msg.str());
  }
}

/// \brief Parse "param_chem_pot" and store as
///     "param_chem_pot" vector values
///
/// If successfully parsed, `parser->value` will contain a
/// monte::ValueMap with:
/// - vector_values["param_chem_pot"]: (size = system components size)
///
/// If unsuccesfully parsed, `parser.valid() == false`.
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
/// Requires:
/// - get_composition_converter(system_type const &system);
void parse_param_chem_pot(InputParser<monte::ValueMap> &parser,
                          std::shared_ptr<system_type> const &system) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_param_chem_pot: parser must have non-empty value");
  }

  if (!parser.self.contains("param_chem_pot")) {
    return;
  }

  auto const &composition_converter = get_composition_converter(*system);
  if (parser.self["param_chem_pot"].is_array()) {
    Eigen::VectorXd &value = parser.value->vector_values["param_chem_pot"];
    parser.optional(value, "param_chem_pot");
    if (value.size() != composition_converter.independent_compositions()) {
      std::stringstream msg;
      msg << "Error: 'param_chem_pot' size mismatch.";
      parser.insert_error("param_chem_pot", msg.str());
    }
  } else if (parser.self["param_chem_pot"].is_object()) {
    std::map<std::string, double> input;
    parser.optional(input, "param_chem_pot");
    try {
      parser.value->vector_values["param_chem_pot"] =
          make_param_chem_pot(composition_converter, input);
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
}

/// \brief Parse "include_formation_energy" or set default value
void parse_include_formation_energy(InputParser<monte::ValueMap> &parser) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_include_formation_energy: parser must have non-empty "
        "value");
  }

  bool default_value = true;
  if (parser.self.contains("corr_matching_pot")) {
    default_value = false;
  } else if (parser.self.contains("random_alloy_corr_matching_pot")) {
    default_value = false;
  }

  parser.optional_else<bool>(
      parser.value->boolean_values["include_formation_energy"],
      "include_formation_energy", default_value);
}

/// \brief Parse boolean value to monte::ValueMap
void parse_boolean(InputParser<monte::ValueMap> &parser, std::string option) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_boolean: parser must have non-empty value");
  }

  if (!parser.self.contains(option)) {
    return;
  }
  parser.optional(parser.value->boolean_values[option], option);
}

/// \brief Parse scalar value to monte::ValueMap
void parse_scalar(InputParser<monte::ValueMap> &parser, std::string option) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_scalar: parser must have non-empty value");
  }

  if (!parser.self.contains(option)) {
    return;
  }
  parser.optional(parser.value->scalar_values[option], option);
}

/// \brief Parse vector value to monte::ValueMap
void parse_vector(InputParser<monte::ValueMap> &parser, std::string option) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_vector: parser must have non-empty value");
  }

  if (!parser.self.contains(option)) {
    return;
  }
  parser.optional(parser.value->vector_values[option], option);
}

/// \brief Parse matrix value to monte::ValueMap
void parse_matrix(InputParser<monte::ValueMap> &parser, std::string option) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_matrix: parser must have non-empty value");
  }

  if (!parser.self.contains(option)) {
    return;
  }
  parser.optional(parser.value->matrix_values[option], option);
}

/// \brief Parse "corr_matching_pot"
///
/// Expected format:
/// \code
///   "corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "targets": [
///       {"index":int, "value":number (default=0.), "weight":number
///       (default=1.)},
///       ...
///     ]
///   }
/// \endcode
///
void parse_corr_matching_pot(InputParser<monte::ValueMap> &parser,
                             bool is_increment) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_corr_matching_pot: parser must have non-empty value");
  }

  if (!parser.self.contains("corr_matching_pot")) {
    return;
  }

  try {
    CorrMatchingParams params;
    from_json(params, parser.self["corr_matching_pot"]);
    if (is_increment) {
      parser.value->vector_values["corr_matching_pot"] = to_VectorXd(params);
    } else {
      parser.value->vector_values["corr_matching_pot"] =
          to_VectorXd_increment(params);
    }
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error: could not parse 'corr_matching_pot'.";
    parser.insert_error("corr_matching_pot", e.what());
  }
}

/// \brief Parse "random_alloy_corr_matching_pot"
///
/// Expected format:
/// \code
///   "random_alloy_corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "sublattice_prob": [
///       [number, number, ...], // sublattice 0
///       [number, number, ...], // sublattice 1
///       ...
///     ],
///     "targets": [
///       {"index":int, "weight":number (default=1.)},
///       ...
///     ]
///   }
/// \endcode
///
void parse_random_alloy_corr_matching_pot(
    InputParser<monte::ValueMap> &parser,
    std::shared_ptr<system_type> const &system, bool is_increment) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error in parse_random_alloy_corr_matching_pot: parser must have "
        "non-empty value");
  }

  if (!parser.self.contains("random_alloy_corr_matching_pot")) {
    return;
  }

  try {
    RandomAlloyCorrMatchingParams params(get_random_alloy_corr_f(*system));
    from_json(params, parser.self["random_alloy_corr_matching_pot"]);
    if (is_increment) {
      parser.value->vector_values["random_alloy_corr_matching_pot"] =
          to_VectorXd(params);
    } else {
      parser.value->vector_values["random_alloy_corr_matching_pot"] =
          to_VectorXd_increment(params);
    }
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "Error: could not parse 'random_alloy_corr_matching_pot'.";
    parser.insert_error("random_alloy_corr_matching_pot", e.what());
  }
}

}  // namespace clexmonte
}  // namespace CASM
