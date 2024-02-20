#ifndef CASM_clexmonte_FixedConfigGenerator
#define CASM_clexmonte_FixedConfigGenerator

#include <vector>

#include "casm/clexmonte/definitions.hh"
#include "casm/clexmonte/run/ConfigGenerator.hh"
#include "casm/clexmonte/run/RunData.hh"
#include "casm/monte/run_management/State.hh"

namespace CASM {
namespace clexmonte {

/// \brief A `ConfigGenerator` for state generation -- always returns the same
/// configuration
///
/// - Returns the same configuration no matter what the current
///   conditions and previous runs are.
class FixedConfigGenerator : public ConfigGenerator {
 public:
  FixedConfigGenerator(config_type const &configuration)
      : m_configuration(configuration) {}

  config_type operator()(monte::ValueMap const &conditions,
                         std::vector<RunData> const &completed_runs) override {
    return m_configuration;
  }

 private:
  config_type m_configuration;
};

}  // namespace clexmonte
}  // namespace CASM

#endif
