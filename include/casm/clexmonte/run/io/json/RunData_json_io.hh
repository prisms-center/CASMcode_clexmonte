#ifndef CASM_clexmonte_RunData_json_io
#define CASM_clexmonte_RunData_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clexmonte/run/RunData.hh"
#include "casm/clexmonte/state/io/json/State_json_io.hh"
#include "casm/monte/io/json/ValueMap_json_io.hh"

namespace CASM {

inline jsonParser &to_json(clexmonte::RunData const &run_data, jsonParser &json,
                           bool write_initial_states, bool write_final_states) {
  if (write_initial_states) {
    json["initial_state"] = run_data.initial_state;
  }
  if (write_final_states) {
    json["final_state"] = run_data.final_state;
  }
  json["conditions"] = run_data.conditions;
  json["transformation_matrix_to_supercell"] =
      run_data.transformation_matrix_to_super;
  json["n_unitcells"] = run_data.n_unitcells;
  return json;
}

inline void parse(InputParser<clexmonte::RunData> &parser,
                  config::SupercellSet &supercells) {
  parser.value = std::make_unique<clexmonte::RunData>();

  clexmonte::RunData &run_data = *parser.value;
  parser.optional(run_data.initial_state, "initial_state", supercells);
  parser.optional(run_data.final_state, "final_state", supercells);
  parser.require(run_data.conditions, "conditions");
  parser.require(run_data.transformation_matrix_to_super,
                 "transformation_matrix_to_supercell");
  parser.require(run_data.n_unitcells, "n_unitcells");
}

inline void from_json(clexmonte::RunData &run_data, jsonParser const &json,
                      config::SupercellSet &supercells) {
  InputParser<clexmonte::RunData> parser{json, supercells};

  std::runtime_error error_if_invalid{
      "Error reading clexmonte::RunData from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  run_data = *parser.value;
}

inline jsonParser &to_json(clexmonte::RunDataOutputParams &output_params,
                           jsonParser &json) {
  json["save_all_initial_states"] = output_params.do_save_all_initial_states;
  json["save_all_final_states"] = output_params.do_save_all_final_states;
  json["save_last_final_state"] = output_params.do_save_last_final_state;
  json["write_initial_states"] = output_params.write_initial_states;
  json["write_final_states"] = output_params.write_final_states;
  if (!output_params.output_dir.empty()) {
    json["output_dir"] = output_params.output_dir.string();
  }
  return json;
}

inline void from_json(clexmonte::RunDataOutputParams &output_params,
                      jsonParser const &json) {
  InputParser<clexmonte::RunDataOutputParams> parser{json};

  std::runtime_error error_if_invalid{
      "Error reading clexmonte::RunDataOutputParams from JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  output_params = *parser.value;
}

inline void parse(InputParser<clexmonte::RunDataOutputParams> &parser) {
  parser.value = std::make_unique<clexmonte::RunDataOutputParams>();

  clexmonte::RunDataOutputParams &output_params = *parser.value;
  parser.optional_else(output_params.do_save_all_initial_states,
                       "save_all_initial_states", false);
  parser.optional_else(output_params.do_save_all_final_states,
                       "save_all_final_states", false);
  parser.optional_else(output_params.do_save_last_final_state,
                       "save_last_final_state", true);
  parser.optional_else(output_params.write_initial_states,
                       "write_initial_states", false);
  parser.optional_else(output_params.write_final_states, "write_final_states",
                       false);
  std::string output_dir;
  parser.optional(output_dir, "output_dir");
  output_params.output_dir = fs::path(output_dir);
}

}  // namespace CASM

#endif
