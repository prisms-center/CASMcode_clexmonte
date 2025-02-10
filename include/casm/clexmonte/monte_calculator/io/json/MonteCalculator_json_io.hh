#ifndef CASM_clexmonte_MonteCalculator_json_io
#define CASM_clexmonte_MonteCalculator_json_io

#include <memory>
#include <vector>

#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"
#include "casm/global/filesystem.hh"

namespace CASM {
template <typename T>
class InputParser;
class jsonParser;

namespace clexmonte {
class System;

/// \brief Parse MonteCalculator from JSON
void parse(InputParser<std::shared_ptr<MonteCalculator>> &parser,
           std::shared_ptr<System> &system, jsonParser const &params,
           std::shared_ptr<MonteCalculator::engine_type> engine,
           std::vector<fs::path> search_path = {});

}  // namespace clexmonte
}  // namespace CASM

#endif
