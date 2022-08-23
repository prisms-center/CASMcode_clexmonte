#ifndef CASM_clexmonte_clex_ClexData_json_io
#define CASM_clexmonte_clex_ClexData_json_io

#include <memory>

namespace CASM {

template <typename T>
class InputParser;

namespace clexulator {
class PrimNeighborList;
}

namespace clexmonte {
struct ClexData;
struct LocalClexData;

void parse(InputParser<clexmonte::ClexData> &parser,
           std::shared_ptr<clexulator::PrimNeighborList> &prim_neighbor_list);

void parse(InputParser<clexmonte::LocalClexData> &parser,
           std::shared_ptr<clexulator::PrimNeighborList> &prim_neighbor_list);

}  // namespace clexmonte
}  // namespace CASM

#endif
