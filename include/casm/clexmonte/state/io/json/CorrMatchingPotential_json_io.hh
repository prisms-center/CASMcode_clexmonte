#ifndef CASM_clexmonte_state_CorrMatchingPotential_json_io
#define CASM_clexmonte_state_CorrMatchingPotential_json_io

namespace CASM {
class jsonParser;

namespace clexmonte {

struct CorrMatchingTarget;
struct CorrMatchingParams;
struct RandomAlloyCorrMatchingParams;

jsonParser &to_json(CorrMatchingTarget const &target, jsonParser &json);

void from_json(CorrMatchingTarget &target, jsonParser const &json);

jsonParser &to_json(CorrMatchingParams const &params, jsonParser &json);

void from_json(CorrMatchingParams &params, jsonParser const &json);

jsonParser &to_json(RandomAlloyCorrMatchingParams const &params,
                    jsonParser &json);

void from_json(RandomAlloyCorrMatchingParams &params, jsonParser const &json);

}  // namespace clexmonte
}  // namespace CASM

#endif
