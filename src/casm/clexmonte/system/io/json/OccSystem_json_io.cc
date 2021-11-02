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
  std::cout << "begin parse OccSystem" << std::endl;

  std::cout << "parse OccSystem 1" << std::endl;
  // 1) Read prim
  std::shared_ptr<xtal::BasicStructure const> shared_prim =
      parser.require<xtal::BasicStructure>("prim", TOL);

  std::cout << "parse OccSystem 2" << std::endl;
  // 2) Read composition axes
  std::unique_ptr<composition::CompositionConverter> composition_axes =
      parser.require<composition::CompositionConverter>("composition_axes");

  std::cout << "parse OccSystem 3" << std::endl;
  // 3) Construct formation energy clexulator and coefficients
  std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list;
  auto formation_energy_clex_data_subparser =
      parser.subparse<ClexData>("formation_energy", prim_neighbor_list);

  std::cout << "parse OccSystem 4" << std::endl;
  if (parser.valid()) {
    parser.value = std::make_unique<OccSystem>(
        shared_prim, *composition_axes,
        *formation_energy_clex_data_subparser->value);
  }
  std::cout << "finish parse OccSystem" << std::endl;
}

}  // namespace clexmonte
}  // namespace CASM
