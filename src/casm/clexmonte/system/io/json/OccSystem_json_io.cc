#include "casm/clexmonte/system/io/json/OccSystem_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/io/json/ClexData_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse equivalents_info.json
///
/// TODO: format
void parse(
      InputParser<EquivalentsInfo> &parser,
      config::Prim const &prim) {

  auto const &basicstructure = *prim.basicstructure;

  parser.value = std::make_unique<EquivalentsInfo>();
  auto &info = *parser.value;

  parser.require(info.equivalent_generating_op_indices, "equivalent_generating_ops");

  if (parser.self.contains("equivalents")) {
    auto begin = parser.self["equivalents"].begin();
    auto end = parser.self["equivalents"].end();
    for (auto it=begin; it != end; ++it) {
      auto subparser = parser.subparse<clust::IntegralCluster>("phenomenal", basicstructure);
      if (subparser.valid()) {
        info.phenomenal_clusters.push_back(*subparse->value);
      }
    }
  }

  if (info.equivalent_generating_op_indices.size() != info.phenomenal_clusters.size()) {
    parser.insert_error("equivalent_generating_ops", "Size mismatch with 'equivalents'");
  }

  if (info.equivalent_generating_op_indices.size() == 0) {
    parser.insert_error("equivalent_generating_ops", "Size==0");
  }

  if (!parser.valid()) {
    return;
  }

  for (Index i=0; i<info.equivalent_generating_op_indices.size(); ++i) {
    Index fg_index = info.equivalent_generating_op_indices[i];
    xtal::SymOp equivalence_map_op = make_equivalence_map_op(
        info.phenomenal_clusters[0],
        info.phenomenal_clusters[i],
        prim.lattice().lat_column_mat(),
        prim.sym_info.factor_group->element[fg_index],
        prim.sym_info.unitcellcoord_symgroup_rep[fg_index]);
    info.equivalent_generating_ops.push_back(equivalence_map_op);
  }
}

/// \brief Parse "events"/<event_name> from JSON
void parse(
      InputParser<OccEventTypeDescription> &parser,
      config::Prim const &prim,
      std::map<std::string, EquivalentsInfo> const &equivalents_info,
      occ_events::OccSystem const &event_system) {

  // parse "kra"
  std::string kra_clex_name;
  parser.require(kra_clex_name, "kra");
  if (!equivalents_info.count(kra_clex_name)) {
    parser.insert_error(fs::path("events") / it.name() / "kra", "KRA local clex not found")
  }

  // parse "freq"
  std::string freq_clex_name;
  parser.require(freq_clex_name, "freq");
  if (!equivalents_info.count(freq_clex_name)) {
    parser.insert_error(fs::path("events") / it.name() / "freq", "Attempt frequency local clex not found")
  }

  // parse "event"
  std::string _filename;
  parser.require(_filename, "event");

  if (!parser.valid()) {
    return;
  }

  jsonParser json(fs::path(_filename));
  InputParser<OccEvent> parser(json, event_system);

  if (!parser.valid()) {
    return;
  }

  parser.value = std::make_unique<OccEventTypeDescription>();
  event_type = *parser.value;
  event_type.kra_clex_name = kra_clex_name;
  event_type.freq_clex_name = freq_clex_name;

  // generate equivalent events, consistent with kra clex
  // (freq clex must also be consistent)
  occ_events::OccEvent event = *parser.value;
  auto const &info = equivalents_info.at(kra_clex_name);
  std::vector<OccEventRep> occevent_reps = make_occevent_symgroup_rep(
    info.equivalent_generating_ops, basicstructure);
  for (auto const &rep  : occevent_reps) {
    event_type.events.push_back(copy_apply(rep, event));
  }
}


/// \brief Parse OccSystem from JSON
void parse(InputParser<OccSystem> &parser) {

  // Parse prim
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      parser.require<xtal::BasicStructure>("prim", TOL);

  // Parse composition axes
  std::unique_ptr<composition::CompositionConverter> composition_axes =
      parser.require<composition::CompositionConverter>("composition_axes");

  if (!parser.valid()) {
    return;
  }

  // Construct OccSystem
  parser.value = std::make_unique<OccSystem>(shared_prim, *composition_axes);
  auto &system = *parser.value;
  auto const &prim = *system.prim;
  auto &prim_neighbor_list = system.prim_neighbor_list;

  // Parse basis_sets
  if (parser.self.contains("basis_sets")) {
    auto &basis_sets = system.basis_sets;
    auto begin = parser.self["basis_sets"].begin();
    auto end = parser.self["basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      std::string source;
      parser.require(source, fs::path("basis_sets") / it.name() / "source");
      auto subparser = parser.subparse<clexulator::Clexulator>(
          fs::path("basis_sets") / it.name(), prim_neighbor_list);
      if (subparser.valid()) {
        auto clexulator = std::make_shared<clexulator::Clexulator>(
            std::move(*subparser->value));
        basis_sets.emplace(source, clexulator);
      }
    }
  }

  // Parse local_basis_sets && equivalents_info
  if (parser.self.contains("local_basis_sets")) {

    // parse to construct local Clexulator
    auto &local_basis_sets = system.local_basis_sets;
    auto begin = parser.self["local_basis_sets"].begin();
    auto end = parser.self["local_basis_sets"].end();
    for (auto it = begin; it != end; ++it) {
      auto subparser = parser.subparse<std::vector<clexulator::Clexulator>>(
          fs::path("local_basis_sets") / it.name(), prim_neighbor_list);
      if (subparser.valid()) {
        auto local_clexulator =
            std::make_shared<std::vector<clexulator::Clexulator>>(
                std::move(*subparser->value));
        basis_sets.emplace(it.name(), local_clexulator);
      }
    }

    // parse equivalents_info files
    auto &equivalents_info = system.equivalents_info;
    auto begin = parser.self["local_basis_sets"].begin();
    auto end = parser.self["local_basis_sets"].end();
    for (auto it = begin; it != end; ++it) {

      std::string _filename;
      parser.require(_filename, fs::path("local_basis_sets") / it.name() / "equivalents_info");

      if (!fs::exists(_filename)) {
        parser.insert_error(opt, "No equivalents_info.json file found");
        continue;
      }

      jsonParser json(fs::path(_filename));
      InputParser<EquivalentsInfo> info_parser(json, *prim);
      if (info_parser.valid()) {
        equivalents_info.emplace(it.name(), std::move(*info_parser.value));
      }
    }
  }

  // Parse clex
  if (parser.self.contains("clex")) {
    auto &clex_data = system.clex_data;
    auto begin = parser.self["clex"].begin();
    auto end = parser.self["clex"].end();
    for (auto it = begin; it != end; ++it) {
      std::string basis_set_name;
      parser.require(basis_set_name,
                     fs::path("clex") / it.name() / "basis_set");
      auto basis_set_it = system.basis_set.find(basis_set_name);
      if (basis_set_it == system.basis_set.end()) {
        parser.insert_error(fs::path("clex") / it.name(),
                            "No matching basis_set");
        continue;
      }
      auto subparser =
          parser.subparse<ClexData>(fs::path("clex") / it.name(),
                                    prim_neighbor_list, basis_set_it->second);
      if (subparser.valid()) {
        clex_data.emplace(it.name(), std::move(*subparser->value));
      }
    }
  }

  // Parse local_clex
  if (parser.self.contains("local_clex")) {
    auto &local_clex_data = system.local_clex_data;
    auto begin = parser.self["local_clex"].begin();
    auto end = parser.self["local_clex"].end();
    for (auto it = begin; it != end; ++it) {
      std::string local_basis_set_name;
      parser.require(local_basis_set_name,
                     fs::path("local_clex") / it.name() / "local_basis_set");
      auto local_basis_set_it =
          system.local_basis_set.find(local_basis_set_name);
      if (local_basis_set_it == system.local_basis_set.end()) {
        parser.insert_error(fs::path("local_clex") / it.name(),
                            "No matching local_basis_set");
        continue;
      }
      auto subparser = parser.subparse<LocalClexData>(
          fs::path("clex") / it.name(), prim_neighbor_list,
          local_basis_set_it->second);
      if (subparser.valid()) {
        local_clex_data.emplace(it.name(), std::move(*subparser->value));
      }
    }
  }

  // Parse event_system && events
  if (parser.self.contains("events")) {

    // parse "event_system"
    auto event_system_subparser = parser.subparse<occ_events::OccSystem>("event_system", shared_prim);
    if (event_system_subparser->valid()) {
      system.event_system = std::make_shared<occ_events::OccSystem>(
          std::move(*event_system_subparser.value));
    }

    // parse "events"
    auto &event_data = system.event_types;
    auto begin = parser.self["events"].begin();
    auto end = parser.self["events"].end();
    for (auto it = begin; it != end; ++it) {
      auto subparser = parser.subparse<OccEventTypeDescription>(
          fs::path("events") / it.name(),
          *system.prim, system.equivalents_info, *system.event_system);
      if (subparser.valid()) {
        event_types.emplace(it.name(), std::move(*subparser->value));
      }
    }
  }

  // Parse DoFSpaces
  // TODO...

}

}  // namespace clexmonte
}  // namespace CASM
