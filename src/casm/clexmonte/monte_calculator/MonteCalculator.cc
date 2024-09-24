
#include "casm/clexmonte/monte_calculator/MonteCalculator.hh"

#include <filesystem>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {
namespace clexmonte {

namespace MonteCalculator_impl {

void print_runtime_lib_options_help(std::ostream &sout) {
  sout << "Error compiling MonteCalculator. To fix: \n";
  sout << "  - Check compiler error messages.\n";
  sout << "  - Check compiler options with 'casm settings -l'\n";
  sout << "    - Update compiler options with 'casm settings "
          "--set-compile-options '...options...'\n";
  sout << "    - Make sure the casm headers can be found by including "
          "'-I/path/to/casm'\n";
};

/// Make shared_ptr<RuntimeLibrary>, logging progress and errors
std::shared_ptr<RuntimeLibrary> make_shared_runtime_lib(
    std::string filename_base, std::string compile_options,
    std::string so_options, std::string compile_msg) {
  auto &log = CASM::log();
  bool write_status_messages = false;
  if (!fs::exists(filename_base + ".so")) {
    write_status_messages = true;
  }
  if (write_status_messages) {
    log.compiling<Log::standard>(filename_base + ".cc");
    log.begin_lap();
    log << compile_msg << std::endl;
  }
  try {
    std::shared_ptr<RuntimeLibrary> result = std::make_shared<RuntimeLibrary>(
        filename_base, compile_options, so_options);
    if (write_status_messages) {
      log << "compile time: " << log.lap_time() << " (s)\n" << std::endl;
    }
    return result;
  } catch (runtime_lib_compile_error &e) {
    e.print(err_log());
    print_runtime_lib_options_help(log);
    throw;
  } catch (runtime_lib_shared_error &e) {
    e.print(err_log());
    print_runtime_lib_options_help(log);
    throw;
  } catch (std::exception &e) {
    print_runtime_lib_options_help(log);
    throw;
  }
}

}  // namespace MonteCalculator_impl

/// \brief MonteCalculator factory function
///
/// This does the following:
/// 1. Construct a std::shared_ptr<MonteCalculator> with the BaseMonteCalculator
/// pointer to the implementation instance, and with the runtime library it came
/// from (if applicable).
/// 2. Pass parameters and system data. This triggers checks for required system
/// data and required and optional parameters (existence only). A warning
/// is printed to CASM::log for additional parameters (top level JSON
/// attributes only, excluding those starting with '_').
/// 3. Standard sampling, analysis, and modifying functions are constructed as
/// lambdas which are provided the std::shared_ptr<MonteCalculator> and then
/// added to the calculator. Additional functions can be added as desired
/// before running the method.
///
/// \param params Calculation method parameters, as specified by the
///     particular calculation type
/// \param system System data
/// \param base_calculator The underlying calculator implementation instance
/// \param lib If the `base_calculator` is from a runtime library, it should be
///     provided. Otherwise, give nullptr.
/// \return The std::shared_ptr<MonteCalculator>
std::shared_ptr<MonteCalculator> make_monte_calculator(
    jsonParser const &params, std::shared_ptr<system_type> system,
    std::unique_ptr<BaseMonteCalculator> base_calculator,
    std::shared_ptr<RuntimeLibrary> lib) {
  std::shared_ptr<MonteCalculator> calculator =
      std::make_shared<MonteCalculator>(std::move(base_calculator), lib);

  /// Pass parameters and system data.
  calculator->reset(params, system);

  /// Standard sampling, analysis, and modifying functions are provided the
  /// std::shared_ptr<MonteCalculator> and added to the calculator
  calculator->sampling_functions =
      calculator->standard_sampling_functions(calculator);
  calculator->json_sampling_functions =
      calculator->standard_json_sampling_functions(calculator);
  calculator->analysis_functions =
      calculator->standard_analysis_functions(calculator);
  calculator->modifying_functions =
      calculator->standard_modifying_functions(calculator);

  std::optional<monte::SelectedEventDataFunctions> x =
      calculator->standard_selected_event_data_functions(calculator);
  if (x.has_value()) {
    calculator->selected_event_data_functions = std::move(x.value());
  }

  return calculator;
}

/// \brief MonteCalculator factory function, from source
std::shared_ptr<MonteCalculator> make_monte_calculator_from_source(
    fs::path dirpath, std::string calculator_name, jsonParser const &params,
    std::shared_ptr<system_type> system, std::string compile_options,
    std::string so_options) {
  // Construct the RuntimeLibrary
  std::shared_ptr<RuntimeLibrary> lib;
  try {
    lib = MonteCalculator_impl::make_shared_runtime_lib(
        (dirpath / calculator_name).string(), compile_options, so_options,
        "compile time depends on calculator complexity");
  } catch (std::exception &e) {
    CASM::log()
        << "MonteCalculator construction failed: could not construct runtime "
           "library."
        << std::endl;
    throw;
  }

  // Get the factory function
  std::function<clexmonte::BaseMonteCalculator *(void)> factory;
  factory = lib->get_function<clexmonte::BaseMonteCalculator *(void)>(
      "make_" + calculator_name);

  // Use the factory to construct the BaseMonteCalculator
  std::unique_ptr<clexmonte::BaseMonteCalculator> base(factory());

  // Then use `make_monte_calculator` to construct a shared MonteCalculator, set
  // parameters, and add standard sampling functions
  return make_monte_calculator(params, system, std::move(base), lib);
}

Eigen::VectorXd scalar_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  monte::ValueMap const &conditions =
      calculation->state_data()->state->conditions;
  if (!conditions.scalar_values.count(key)) {
    std::stringstream msg;
    msg << "Error accessing MonteCalculator scalar condition '" << key
        << "': not found";
    throw std::runtime_error(msg.str());
  }
  return monte::reshaped(conditions.scalar_values.at(key));
}

Eigen::VectorXd vector_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  monte::ValueMap const &conditions =
      calculation->state_data()->state->conditions;
  if (!conditions.vector_values.count(key)) {
    std::stringstream msg;
    msg << "Error accessing MonteCalculator vector condition '" << key
        << "': not found";
    throw std::runtime_error(msg.str());
  }
  return conditions.vector_values.at(key);
}

Eigen::VectorXd matrix_conditions(
    std::shared_ptr<MonteCalculator> const &calculation, std::string key) {
  monte::ValueMap const &conditions =
      calculation->state_data()->state->conditions;
  if (!conditions.matrix_values.count(key)) {
    std::stringstream msg;
    msg << "Error accessing MonteCalculator matrix condition '" << key
        << "': not found";
    throw std::runtime_error(msg.str());
  }
  return conditions.matrix_values.at(key);
}

system_type const &get_system(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return *calculation->system();
}

state_type const &get_state(
    std::shared_ptr<MonteCalculator> const &calculation) {
  return *calculation->state_data()->state;
}

/// \brief Make temporary monte::OccLocation if necessary
///
/// \param occ_location Reference-to-pointer. Use the pointed to OccLocation if
///     not nullptr. If nullptr, construct a temporary monte::OccLocation and
///     set `occ_location` to point at it.
/// \param tmp Where to construct temporary monte::OccLocation if
///     `occ_location` is nullptr.
/// \param calculation Where to get data needed to
///     construct temporary monte::OccLocation
void make_temporary_if_necessary(state_type const &state,
                                 monte::OccLocation *&occ_location,
                                 std::unique_ptr<monte::OccLocation> &tmp,
                                 MonteCalculator const &calculation) {
  if (!occ_location) {
    auto const &system_ptr = calculation.system();
    if (!system_ptr) {
      throw std::runtime_error(
          "Error checking if a temporary OccLocation is necessary: "
          "occ_location is null and system is null");
    }
    auto &system = *system_ptr;
    make_temporary_if_necessary(state, occ_location, tmp, system,
                                calculation.update_atoms(),
                                calculation.save_atom_info());
  }
}

}  // namespace clexmonte
}  // namespace CASM
