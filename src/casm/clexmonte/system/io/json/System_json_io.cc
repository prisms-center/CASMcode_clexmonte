#include "casm/clexmonte/system/io/json/System_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/System.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/clexulator/io/json/Clexulator_json_io.hh"
#include "casm/clexulator/io/json/SparseCoefficients_json_io.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

namespace CASM {
namespace clexmonte {

namespace {

/// \brief Parse
template <typename ParserType, typename RequiredType>
bool parse_from_files_object(ParserType &parser, fs::path option,
                             std::vector<RequiredType> &vec,
                             std::map<std::string, Index> &glossary) {
  auto obj_it = parser.self.find(option);
  if (obj_it == parser.self.end() || !obj_it->is_obj()) {
    parser.insert_error(option, "Missing required JSON object");
    return false;
  }
  Index i = 0;
  for (auto it = obj_it->begin(); it != obj_it->end(); ++it) {
    auto subparser = parser.template subparse_from_file<RequiredType>(
        option / std::to_string(i));
    if (!subparser->valid()) {
      return false;
    }
    vec.push_back(std::move(*subparser->value));
    glossary.emplace(it.name(), i);
    ++i;
  }
  return true;
}

/// \brief Parse from file
template <typename ParserType, typename RequiredType, typename... Args>
bool parse_from_file(ParserType &parser, fs::path option, RequiredType &value,
                     Args &&...args) {
  auto subparser = parser.template subparse_from_file<RequiredType>(
      option, std::forward<Args>(args)...);
  if (!subparser->valid()) {
    return false;
  }
  value = std::move(*subparser->value);
  return true;
}

/// \brief Parse and validate "basis_set" or "local_basis_set" (name)
template <typename ParserType, typename BasisSetMapType>
bool parse_and_validate_basis_set_name(ParserType &parser, fs::path option,
                                       std::string &basis_set_name,
                                       BasisSetMapType const &basis_sets) {
  parser.require(basis_set_name, option);

  // validate basis_set_name exists
  auto basis_set_it = basis_sets.find(basis_set_name);
  if (basis_set_it == basis_sets.end()) {
    parser.insert_error(option, "No basis set with matching name");
    return false;
  }
  return true;
}

/// \brief Parse "events"/<event_name> from JSON
///
/// TODO: document format (see
/// tests/unit/clexmonte/data/kmc/system_template.json)
/// - event : occ_events::OccEvent file path
/// - local_basis_set: name (matching one in System::local_basis_sets /
/// equivalents_info)
/// - coefficients/kra: SparseCoefficients file path
/// - coefficients/freq: SparseCoefficients file path
template <typename ParserType>
bool parse_event(
    ParserType &parser, fs::path option,
    std::map<std::string, OccEventTypeData> &event_type_data,
    std::map<std::string, LocalMultiClexData> &local_multiclex,
    xtal::BasicStructure const &prim,
    std::map<std::string,
             std::shared_ptr<std::vector<clexulator::Clexulator>>> const
        &local_basis_sets,
    std::map<std::string, EquivalentsInfo> const &equivalents_info,
    std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep,
    occ_events::OccSystem const &event_system) {
  std::string event_name = option.filename();

  auto coeffs_it = parser.self.find_at(option / "coefficients");
  if (coeffs_it == parser.self.end()) {
    parser.insert_error(option / "coefficients",
                        "Missing coefficients info for this local basis set");
    return false;
  }
  if (!coeffs_it->is_obj()) {
    parser.insert_error(
        option / "coefficients",
        "Coefficients info must be a JSON object of <key>:<file path>");
    return false;
  }
  Index coeffs_size = coeffs_it->size();

  // --- Create a local multi-cluster expansion for the event type ---
  LocalMultiClexData curr_local_multiclex;
  curr_local_multiclex.coefficients.resize(coeffs_size);

  // parse option / "coefficients" / <name> : <coefficients file path>
  Index i = 0;
  for (auto it = coeffs_it->begin(); it != coeffs_it->end(); ++it) {
    std::string key = it.name();
    if (!parse_from_file(parser, option / "coefficients" / key,
                         curr_local_multiclex.coefficients[i])) {
      return false;
    }
    curr_local_multiclex.coefficients_glossary.emplace(key, i);
    ++i;
  }

  // parse and validate "local_basis_set"
  if (!parse_and_validate_basis_set_name(
          parser, option / "local_basis_set",
          curr_local_multiclex.local_basis_set_name, local_basis_sets)) {
    return false;
  }
  if (equivalents_info.find(curr_local_multiclex.local_basis_set_name) ==
      equivalents_info.end()) {
    parser.insert_error(option / "local_basis_set",
                        "Missing equivalents_info for this local basis set");
    return false;
  }

  // Save the LocalMultiClexData
  local_multiclex.emplace(event_name, curr_local_multiclex);

  /// --- Parse "event", construct and validate equivalents ---

  // parse "event" (from file)
  auto event_subparser =
      parser.template subparse_from_file<occ_events::OccEvent>(option / "event",
                                                               event_system);
  if (!event_subparser->valid()) {
    return false;
  }

  // construct event_data using default constructor
  OccEventTypeData curr_event_type_data;
  curr_event_type_data.local_multiclex_name = event_name;

  // generate equivalent events
  occ_events::OccEvent const &event = *event_subparser->value;
  EquivalentsInfo const &info =
      equivalents_info.at(curr_local_multiclex.local_basis_set_name);
  curr_event_type_data.events =
      make_equivalents(event, info, occevent_symgroup_rep);

  // double-check consistency of event phenomenal clusters and
  // the local basis set phenomenal clusters
  if (!is_same_phenomenal_clusters(curr_event_type_data.events, info)) {
    parser.insert_error(
        option / "local_basis_set",
        "Error generating equivalent events. There is a mismatch between the "
        "event and the phenomenal clusters of the local basis set.");
    return false;
  }

  event_type_data.emplace(event_name, curr_event_type_data);
  return true;
}

}  // namespace

/// \brief Parse System from JSON
///
/// Expected format:
/// \code
///   "prim": <xtal::BasicStructure or file path>
///       Specifies the primitive crystal structure and allowed DoF. Must
///       be the prim used to generate the cluster expansion.
///
///   "composition_axes": <composition::CompositionConverter>
///       Specifies composition axes
///
///   "basis_sets": object (optional)
///       Input specifies one or more CASM cluster expansion basis sets. A JSON
///       object containing one or more of:
///
///           "<name>": {
///             "source": "<path>"
///           }
///
///        where "<name>" is the basis set name, and "source" is the path to a
///        CASM clexulator source file, i.e.
///
///            "/path/to/basis_sets/bset.energy/ZrO_Clexulator_energy.cc"
///
///   "local_basis_sets": object (optional)
///       Input specifies one or more CASM local-cluster expansion basis sets.
///.      A JSON object containing one or more of:
///
///           "<name>": {
///             "source": "<path>",
///             "equivalents_info": "<path">
///           }
///
///        where "<name>" is the local-cluster basis set name, and "source" is
///        the path to a CASM clexulator source file, i.e.
///
///            "/path/to/basis_sets/bset.local/ZrO_Clexulator_local.cc",
///
///        and "equivalents_info" is the path to a "equivalents_info.json" file,
///        i.e.
///
///            "/path/to/basis_sets/bset.local/equivalents_info.json"
///
///   "clex": object (optional)
///       Input specifies one or more cluster expansions. A JSON object
///       containing one or more of:
///
///           "<name>": {
///             "basis_set": "<basis set name>",
///             "coefficients": "<path>"
///           }
///
///        where "<name>" is the cluster expansion name, "<basis set name>" is
///        the name of the basis set for the cluster expansion (matching a key
///        in "basis_sets"), and "coefficients" is the path to a CASM
///        cluster expansion coefficients file.
///
///   "multiclex": object (optional)
///       Input specifies one or more cluster expansions which share basis sets.
///       Similar to "clex", but "coefficients" is a JSON object containing
///       one or more coefficients files:
///
///           "<name>": {
///             "basis_set": "<basis set name>",
///             "coefficients": {
///               "<key>": "<path>"
///               "<key>": "<path>",
///               ...
///           }
///
///   "local_clex": object (optional)
///       Input specifies one or more local-cluster expansions. A JSON object
///       containing one or more of:
///
///           "<name>": {
///             "local_basis_set": "<local basis set name>",
///             "coefficients": "<path>"
///           }
///
///        where "<name>" is the cluster expansion name, "<local basis set
///        name>" is the name of the local-cluster basis set for the
///        local-cluster expansion (matching a key in "local_basis_sets"), and
///        "coefficients" is the path to a CASM cluster expansion coefficients
///        file.
///
///   "local_multiclex": object (optional)
///       Input specifies one or more local-cluster expansions which share
///       basis sets. Similar to "local_clex", but "coefficients" is a JSON
///       object containing one or more coefficients files:
///
///           "<name>": {
///             "local_basis_set": "<local basis set name>",
///             "coefficients": {
///               "<key>": "<path>"
///               "<key>": "<path>",
///               ...
///           }
///
///   "events": object (optional)
///       Input specifies KMC events. A JSON object specifiying one or
///       more events:
///
///           "<event name>": {
///             "event": "<event description file path>",
///             "local_basis_set": "<local basis set name>",
///             "coefficients": {
///               "kra": "<kra path>",
///               "freq": "<freq path>"
///             }
///           }
///
///       where "<event name>" is a name given to the event,
///       "<event description file path>" is the path to a file defining a KMC
///       event (CASM::occ_events::OccEvent in JSON format), "<local basis set
///       name>" is the name of the local-cluster basis set for the
///       local-cluster expansion (matching a key in "local_basis_sets"), "<kra
///       path>" is the path to a cluster expansion coefficients file for the
///       event KRA, and
///       "<freq path>" is the path to a cluster expansion coefficients file
///       for the event attempt frequency. A local multi-clex with name matching
///       "<event name>" will be generated.
///
///   "event_system": string (optional)
///       Input specifies the path to a file specifying how some KMC events are
///       defined (a CASM::occ_events::OccSystem in JSON format).
///
/// \endcode
///
void parse(InputParser<System> &parser) {
  // Parse "prim"
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      parser.require<xtal::BasicStructure>("prim", TOL);

  // Parse "composition_axes"
  std::unique_ptr<composition::CompositionConverter> composition_axes =
      parser.require<composition::CompositionConverter>("composition_axes");

  if (!parser.valid()) {
    return;
  }

  // Construct System
  parser.value = std::make_unique<System>(shared_prim, *composition_axes);
  System &system = *parser.value;

  // Parse "basis_sets"
  if (parser.self.contains("basis_sets")) {
    auto &basis_sets = system.basis_sets;
    auto &prim_neighbor_list = system.prim_neighbor_list;

    auto begin = parser.self["basis_sets"].begin();
    auto end = parser.self["basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      // parse "basis_sets"/<name>/"source"
      auto subparser = parser.subparse<clexulator::Clexulator>(
          fs::path("basis_sets") / it.name(), prim_neighbor_list);
      if (subparser->valid()) {
        auto clexulator = std::make_shared<clexulator::Clexulator>(
            std::move(*subparser->value));
        basis_sets.emplace(it.name(), clexulator);
      }
    }
  }

  // Parse "local_basis_sets"
  if (parser.self.contains("local_basis_sets")) {
    auto &local_basis_sets = system.local_basis_sets;
    auto &equivalents_info = system.equivalents_info;
    auto &prim_neighbor_list = system.prim_neighbor_list;
    auto &prim = *system.prim;

    // construct local Clexulator
    auto begin = parser.self["local_basis_sets"].begin();
    auto end = parser.self["local_basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      // parse "local_basis_sets"/<name>/"source"
      auto subparser = parser.subparse<std::vector<clexulator::Clexulator>>(
          fs::path("local_basis_sets") / it.name(), prim_neighbor_list);
      if (subparser->valid()) {
        auto local_clexulator =
            std::make_shared<std::vector<clexulator::Clexulator>>(
                std::move(*subparser->value));
        local_basis_sets.emplace(it.name(), local_clexulator);
      }

      // parse "local_basis_sets"/<name>/"equivalents_info"
      auto info_subparser = parser.subparse_from_file<EquivalentsInfo>(
          fs::path("local_basis_sets") / it.name() / "equivalents_info", prim);
      if (info_subparser->valid()) {
        equivalents_info.emplace(it.name(), std::move(*info_subparser->value));
      }
    }
  }

  // Parse "clex"
  if (parser.self.contains("clex")) {
    auto &prim = *system.prim;
    auto &clex_data = system.clex_data;
    auto const &basis_sets = system.basis_sets;

    auto begin = parser.self["clex"].begin();
    auto end = parser.self["clex"].end();
    for (auto it = begin; it != end; ++it) {
      ClexData curr;
      fs::path clex_path = fs::path("clex") / it.name();

      // "clex"/<name>/"basis_set"
      if (!parse_and_validate_basis_set_name(parser, clex_path / "basis_set",
                                             curr.basis_set_name, basis_sets)) {
        continue;
      }

      // "clex"/<name>/"coefficients" (SparseCoefficients)
      if (!parse_from_file(parser, clex_path / "coefficients",
                           curr.coefficients)) {
        continue;
      }

      // "clex"/<name>/"coefficients" (BasisSetClusterInfo)
      if (!parse_from_file(parser, clex_path / "coefficients",
                           curr.cluster_info, prim, basis_sets)) {
        continue;
      }

      clex_data.emplace(it.name(), curr);
    }
  }

  // Parse "multiclex"
  if (parser.self.contains("multiclex")) {
    auto &prim = *system.prim;
    auto &multiclex_data = system.multiclex_data;
    auto const &basis_sets = system.basis_sets;

    auto begin = parser.self["multiclex"].begin();
    auto end = parser.self["multiclex"].end();
    for (auto it = begin; it != end; ++it) {
      MultiClexData curr;
      fs::path clex_path = fs::path("multiclex") / it.name();

      // "multiclex"/<name>/"basis_set"
      if (!parse_and_validate_basis_set_name(parser, clex_path / "basis_set",
                                             curr.basis_set_name, basis_sets)) {
        continue;
      }

      // "multiclex"/<name>/"coefficients"
      if (!parse_from_files_object(parser, clex_path / "coefficients",
                                   curr.coefficients,
                                   curr.coefficients_glossary)) {
        continue;
      }

      // "multiclex"/<name>/"coefficients" (BasisSetClusterInfo)
      if (!parse_from_file(parser, clex_path / "coefficients",
                           curr.cluster_info, prim, basis_sets)) {
        continue;
      }

      multiclex_data.emplace(it.name(), curr);
    }
  }

  // Parse "local_clex"
  if (parser.self.contains("local_clex")) {
    auto &local_clex_data = system.local_clex_data;
    auto const &local_basis_sets = system.local_basis_sets;

    auto begin = parser.self["local_clex"].begin();
    auto end = parser.self["local_clex"].end();
    for (auto it = begin; it != end; ++it) {
      LocalClexData curr;
      fs::path clex_path = fs::path("local_clex") / it.name();

      // "local_clex"/<name>/"local_basis_set"
      if (!parse_and_validate_basis_set_name(
              parser, clex_path / "local_basis_set", curr.local_basis_set_name,
              local_basis_sets)) {
        continue;
      }

      // "local_clex"/<name>/"coefficients"
      if (!parse_from_file(parser, clex_path / "coefficients",
                           curr.coefficients)) {
        continue;
      }

      local_clex_data.emplace(it.name(), curr);
    }
  }

  // Parse "local_multiclex"
  if (parser.self.contains("local_multiclex")) {
    auto &local_multiclex_data = system.local_multiclex_data;
    auto const &local_basis_sets = system.local_basis_sets;

    auto begin = parser.self["local_multiclex"].begin();
    auto end = parser.self["local_multiclex"].end();
    for (auto it = begin; it != end; ++it) {
      LocalMultiClexData curr;
      fs::path clex_path = fs::path("local_multiclex") / it.name();

      // "local_multiclex"/<name>/"local_basis_set"
      if (!parse_and_validate_basis_set_name(
              parser, clex_path / "local_basis_set", curr.local_basis_set_name,
              local_basis_sets)) {
        continue;
      }

      // "local_multiclex"/<name>/"coefficients"
      if (!parse_from_files_object(parser, clex_path / "coefficients",
                                   curr.coefficients,
                                   curr.coefficients_glossary)) {
        continue;
      }

      local_multiclex_data.emplace(it.name(), curr);
    }
  }

  // Parse "event_system"
  if (parser.self.contains("event_system")) {
    auto const &basicstructure = system.prim->basicstructure;
    auto event_system_subparser =
        parser.subparse_from_file<occ_events::OccSystem>("event_system",
                                                         basicstructure);
    if (event_system_subparser->valid()) {
      system.event_system = std::make_shared<occ_events::OccSystem>(
          std::move(*event_system_subparser->value));
    }
  }

  // Parse "event_system" and "events"
  if (parser.self.contains("events")) {
    if (system.event_system == nullptr) {
      parser.insert_error("events", "event_system is required to parse events");
    } else {
      // parse "events"/<name>
      auto const &basicstructure = system.prim->basicstructure;
      auto begin = parser.self["events"].begin();
      auto end = parser.self["events"].end();
      for (auto it = begin; it != end; ++it) {
        parse_event(parser, fs::path("events") / it.name(),
                    system.event_type_data, system.local_multiclex_data,
                    *basicstructure, system.local_basis_sets,
                    system.equivalents_info, system.occevent_symgroup_rep,
                    *system.event_system);
      }
    }
  }

  // Parse DoFSpaces
  // TODO...
}

}  // namespace clexmonte
}  // namespace CASM
