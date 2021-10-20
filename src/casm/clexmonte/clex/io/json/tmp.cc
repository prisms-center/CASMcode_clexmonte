
// void parse_composition_conditions(
//     ParentInputParser &parser,
//     fs::path option,
//     std::unique_ptr<VectorValueMap> &conditions,
//     composition::CompositionConverter const &composition_converter,
//     bool as_comp_x) {
//   if (conditions == nullptr) {
//     conditions = std::make_unique<VectorValueMap>();
//   }
//
//   auto it = parser.self.find(option);
//   if (it == parser.self.end()) {
//     parser.insert_error(...)
//   }
//   jsonParser &json = *it;
//
//   auto is_invalid_comp_n_size = [&](jsonParser const &json, int comp_n_size)
//   {
//     if (json.size() != comp_n_size) {
//       parser.insert_error(option / "comp_n", "Mismatch between composition
//       axes and composition input dimension"); return true;
//     }
//     return false;
//   };
//   auto is_invalid_comp_x_size = [&](jsonParser const &json, int comp_x_size)
//   {
//     if (json.size() != comp_x_size) {
//       parser.insert_error(option / "comp_x", "Mismatch between composition
//       axes and composition input dimension"); return true;
//     }
//     return false;
//   }
//   auto is_missing_component = [&](
//       std::string component_name,
//       std::map<std::string, double> const &comp_n_map) {
//     if (!comp_n_map.count(component_name)) {
//       std::stringstream msg;
//       msg << "Invalid component name: " << component_name;
//       parser.insert_error(option / "comp_n", msg.str());
//     }
//     return false;
//   }
//
//
//   int comp_x_size = composition_converter.independent_compositions();
//   int comp_n_size = composition_converter.components().size();
//   it = conditions.find("comp_n");
//   if (it != json.end()) {
//
//     // allowed option: "comp_n": [...array...]
//     if (it->is_array()) {
//       if (is_invalid_comp_n_size(*it, comp_n_size)) { return; }
//       Eigen::VectorXd comp_n;
//       parser.optional(comp_n, "comp_n");
//       emplace_comp_n(comp_n)
//     }
//
//     // allowed option: "comp_n": {<name>:value, <name>:value, ...}
//     else if (it->is_obj()) {
//       if (is_invalid_comp_n_size(*it, comp_n_size)) { return; }
//       std::map<std::string, double> comp_n_map;
//       parser.optional(comp_n_map, "comp_n");
//       Eigen::VectorXd comp_n;
//       int i=0;
//       for (auto const &component_name: composition_converter.components()) {
//         if (is_missing_component(component_name, comp_n_map)) { return; }
//         comp_n(i) = comp_n_map.at(component_name);
//         ++i;
//       }
//     }
//
//     // allowed option: "comp_n": "from_initial_configuration"
//     else if (it->is_string()) {
//       if (is_invalid_comp_n_string(*it, comp_n_size)) { return; }
//
//     }
//
//     else {
//
//     }
//   }
//
//   if (conditions.contains("comp_x")) {
//     auto it = conditions.find("comp_x");
//     if (it->is_array()) {
//
//     }
//     else if (it->is_obj()) {
//
//     }
//     else {
//
//     }
//   }
//
//   parser.insert_error(
//       option,
//       "Error reading composition conditions: Does not contain 'comp_x' or
//       'comp_n'");
// }
// /// \brief Parse canonical Monte Carlo onditions from JSON
// ///
// /// comp_n:
// std::unique_ptr<VectorValueMap>
// parse_canonical_conditions(
//     ParentInputParser &parser,
//     fs::path option,
//     composition::CompositionConverter const &composition_converter) {
//
//   auto conditions = std::make_unique<VectorValueMap>();
//
//   Eigen::VectorXd comp_n;
//   if (json.contains("comp_x")) {
//     Eigen::VectorXd comp_x;
//     int Nparam = composition_converter.independent_compositions();
//     if (json["comp_x"].is_array()) {
//       if (json["comp_x"].size() != Nparam) {
//         parser.insert_error(
//             option,
//             "Mismatch between composition axes and \"comp_x\" size.");
//       }
//       else {
//         from_json(comp_x, json["comp_x"]);
//         comp_n = composition_converter.mol_composition(comp_x);
//       }
//     } else if (json["comp_x"].is_obj()) {
//       comp_x.resize(Nparam);
//
//       for (int i = 0; i < Nparam; i++) {
//         comp_x(i) =
//         json["comp_x"][CompositionConverter::comp_var(i)].get<double>();
//       }
//     }
//
//   } else if (json.contains("comp_n")) {
//     if (json["comp_n"].is_array()) {
//       from_json(comp_n, json["comp_n"]);
//     } else {
//       auto components = composition_converter.components();
//       comp_n.resize(components.size());
//
//       for (int i = 0; i < components.size(); i++) {
//         comp_n(i) = json["comp_n"][components[i]].get<double>();
//       }
//     }
//
//     comp_n = comp_n / comp_n.sum() * shared_prim->basis().size();
//   } else {
//     parser.insert_error(
//         option,
//         "Error reading canonical conditions: Does not contain 'comp_x' or
//         'comp_n'");
//   }
//   else {
//     parser.insert_error(option, "Does not contain \"comp_n\" or \"comp_x\"");
//   }
//   conditions->emplace("comp_n", comp_n);
//
//   double temperature;
//   parser.require(temperature, "temperature");
//
//
// }
