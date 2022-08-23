#include "casm/clexmonte/system/io/json/OccSystem_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/io/json/ClexData_json_io.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexulator/NeighborList.hh"
#include "casm/composition/io/json/CompositionConverter_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

namespace CASM {
namespace clexmonte {

/// \brief Parse OccSystem from JSON
void parse(InputParser<OccSystem> &parser) {
  // 1) Read prim
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      parser.require<xtal::BasicStructure>("prim", TOL);

  // 2) Read composition axes
  std::unique_ptr<composition::CompositionConverter> composition_axes =
      parser.require<composition::CompositionConverter>("composition_axes");

  // 3) Construct formation energy clexulator and coefficients
  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;
  auto formation_energy_clex_data_subparser =
      parser.subparse<ClexData>("formation_energy", prim_neighbor_list);

  // 4) Construct local clexulator and coefficients
  std::map<LocalClexKey, LocalClexData> local_clex_data;
  if (parser.self.contains("local_clex")) {
    auto it = parser.self["local_clex"].begin();
    auto end = parser.self["local_clex"].end();
    for (; it != end; ++it) {
      std::string key = it.name();
      fs::path opt = fs::path("local_clex") / key;
      auto local_clex_data_subparser =
          parser.subparse<LocalClexData>(opt, prim_neighbor_list);
      if (local_clex_data_subparser->valid()) {
        local_clex_data.emplace(key, *local_clex_data_subparser->value);
      }
    }
  }

  // 4) Construct DoFSpaces
  // TODO...

  if (parser.valid()) {
    parser.value = std::make_unique<OccSystem>(
        shared_prim, *composition_axes,
        *formation_energy_clex_data_subparser->value, local_clex_data);
  }
}

}  // namespace clexmonte
}  // namespace CASM
