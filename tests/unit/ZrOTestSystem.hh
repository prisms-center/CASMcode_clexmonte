#ifndef CASM_unittest_ZrOTestSystem
#define CASM_unittest_ZrOTestSystem

#include <filesystem>

#include "autotools.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/system/OccSystem.hh"
#include "casm/clexmonte/system/io/json/OccSystem_json_io.hh"
#include "casm/global/filesystem.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "gtest/gtest.h"
#include "testdir.hh"

namespace test {

/// configure with (in order of priority):
/// - CASM_CXX, CXX, default="c++"
/// - CASM_CXXFLAGS, default="-O3 -Wall -fPIC --std=c++17"
/// - CASM_SOFLAGS, default="-shared"
class ZrOTestSystem : public testing::Test {
 protected:
  ZrOTestSystem() {
    fs::path test_data_dir = test::data_dir("clexmonte") / "ZrOTestSystem";

    fs::path clexulator_src = test_data_dir / "basis_sets" /
                              "bset.formation_energy" /
                              "ZrO_Clexulator_formation_energy.cc";
    fs::path coefficients_src = test_data_dir / "formation_energy_eci.json";

    jsonParser system_json(test_data_dir / "system.json");
    system_json["formation_energy"]["source"] = clexulator_src.string();
    system_json["formation_energy"]["coefficients"] = coefficients_src.string();

    InputParser<clexmonte::OccSystem> parser(system_json);
    std::runtime_error error_if_invalid{"Error reading ZrOTestSystem data"};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    system_data =
        std::shared_ptr<clexmonte::OccSystem>(std::move(parser.value));
  }

  std::shared_ptr<clexmonte::OccSystem> system_data;
};

}  // namespace test

#endif
