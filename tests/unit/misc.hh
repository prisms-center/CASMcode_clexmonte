#ifndef CASM_unittest_misc
#define CASM_unittest_misc

#include <filesystem>
#include <string>

#include "casm/global/definitions.hh"
#include "casm/global/filesystem.hh"

using namespace CASM;

namespace test {

inline Index file_count(fs::path dir) {
  Index count = 0;
  for (auto const& dir_entry : std::filesystem::directory_iterator{dir}) {
    (void)dir_entry;
    ++count;
  }
  return count;
}

}  // namespace test

#endif
