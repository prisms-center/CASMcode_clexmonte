#include "casm/clexmonte/clex/io/json/ClexData_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/clex/ClexData.hh"
#include "casm/clexulator/io/json/SparseCoefficients_json_io.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {
namespace clexmonte {

// - std::shared_ptr<clexulator::PrimNeighborList> prim_neighbor_list
// - fs::path clexulator_src
//   - std::string clexulator_name
//   - fs::path clexulator_dirpath
// - std::string clexulator_compile_options
// - std::string clexulator_so_options
// - fs::path eci_filepath

/// \param parser An InputParser, as genrated by
///     `InputParser::subparse<clexmonte::ClexData>` or one of the other
///     `subparse` methods.
/// \param prim_neighbor_list If not empty, will be used to share neighbor lists
///      among multiple cluster expansions. If sharing, may throw if the
///      Clexulator require different parameters.
/// \param clexulator Clexulator used for cluster expansion basis set
/// evaulation.
///
/// Expected JSON format:
///   coefficients: string (required)
///     Path to a file containing basis set coefficients.
///
void parse(
    InputParser<clexmonte::ClexData> &parser,
    std::shared_ptr<clexulator::PrimNeighborList> const &prim_neighbor_list,
    std::shared_ptr<clexulator::Clexulator> &clexulator) {
  // parse "coefficents"
  std::string _coefficents_path;
  parser.require(_coefficents_path, "coefficients");
  fs::path coefficients_path = _coefficents_path;
  if (!fs::exists(coefficients_path)) {
    parser.insert_error("coefficents",
                        "Error: \"coefficents\" file does not exist.");
  }
  if (!parser.valid()) {
    return;
  }

  // read coefficients file
  jsonParser coefficients_json(coefficients_path);
  auto sparse_coefficients_subparser =
      std::make_shared<InputParser<clexulator::SparseCoefficients>>(
          coefficients_json);
  sparse_coefficients_subparser->type_name =
      CASM::type_name<clexulator::SparseCoefficients>();
  parser.insert("coefficents", sparse_coefficients_subparser);

  if (!parser.valid()) {
    return;
  }

  if (parser.valid()) {
    parser.value = std::make_unique<clexmonte::ClexData>(
        prim_neighbor_list, clexulator, *sparse_coefficients_subparser->value);
  }
}

/// \param parser An InputParser, as genrated by
///     `InputParser::subparse<clexmonte::ClexData>` or one of the other
///     `subparse` methods.
/// \param prim_neighbor_list If not empty, will be used to share neighbor lists
///      among multiple cluster expansions. If sharing, may throw if the
///      Clexulator require different parameters.
/// \param local_clexulator Local Clexulator used for local cluster expansion
/// basis
///     set evaulation.
///
/// Expected JSON format:
///   coefficients: string (required)
///     Path to a file containing basis set coefficients.
void parse(
    InputParser<clexmonte::LocalClexData> &parser,
    std::shared_ptr<clexulator::PrimNeighborList> &prim_neighbor_list,
    std::shared_ptr<std::vector<clexulator::Clexulator>> &local_clexulator) {
  // parse "coefficents"
  std::string _coefficents_path;
  parser.require(_coefficents_path, "coefficients");
  fs::path coefficients_path = _coefficents_path;
  if (!fs::exists(coefficients_path)) {
    parser.insert_error("coefficents",
                        "Error: \"coefficents\" file does not exist.");
  }
  if (!parser.valid()) {
    return;
  }

  // read coefficients file
  jsonParser coefficients_json(coefficients_path);
  auto sparse_coefficients_subparser =
      std::make_shared<InputParser<clexulator::SparseCoefficients>>(
          coefficients_json);
  sparse_coefficients_subparser->type_name =
      CASM::type_name<clexulator::SparseCoefficients>();
  parser.insert("coefficents", sparse_coefficients_subparser);

  // TODO: read "equivalents_info"

  if (parser.valid()) {
    parser.value = std::make_unique<clexmonte::LocalClexData>(
        prim_neighbor_list, local_clexulator,
        *sparse_coefficients_subparser->value);
  }
}

}  // namespace clexmonte
}  // namespace CASM
