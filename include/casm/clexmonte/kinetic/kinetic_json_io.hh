#ifndef CASM_clexmonte_kinetic_json_io
#define CASM_clexmonte_kinetic_json_io

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexmonte/events/io/json/event_data_json_io.hh"
#include "casm/clexmonte/kinetic/kinetic.hh"
#include "casm/clexmonte/misc/parse_array.hh"

namespace CASM {
namespace clexmonte {
namespace kinetic {

/// \brief Parse Kinetic "calculation_options"
///
/// \tparam EngineType
/// \param parser
/// \param system
/// \param random_number_engine (Unused)
///
/// Expected format:
/// \code
///   "event_filters": array of <clexmonte::EventFilterGroup>
///       Allows customizing which events are allowed in which unit cells of
///       the Monte Carlo cell. Each element of the array specifies a set of
///       unit cells for which a particular set of events are allowed. When
///       building the event list, the first filter that a particular unit cell
///       is included in gets applied to that unit cell. If no filter includes
///       a particular unit cell all events are allowed. Each element has the
///       format:
///
///     "exclude": array of int
///         Specifies by index events which will not be allowed.
///     "include": array of int
///         Specifies by index events which will be allowed.
///     "unitcell_index": array of int
///         Specifies by linear unitcell index the unit cells for the
///         events are allowed or not allowed.
///     "include_by_default: bool, (optional, default=true)
///         If `true`, the events not listed explicitly in "include" or
///         "exclude" are allowed. If `false`, the events not listed in
///         "include" or "exclude" are not allowed.
///
/// \endcode
///
template <typename EngineType>
void parse(InputParser<Kinetic<EngineType>> &parser,
           std::shared_ptr<system_type> system,
           std::shared_ptr<EngineType> random_number_engine =
               std::shared_ptr<EngineType>()) {
  // "event_filters"
  std::vector<EventFilterGroup> event_filters;
  if (parser.self.contains("event_filters")) {
    auto subparser = parser.subparse_with<std::vector<EventFilterGroup>>(
        parse_array<EventFilterGroup>, "event_filters", convert);
    if (subparser->valid()) {
      event_filters = std::move(*subparser.value);
    }
  }

  if (parser.valid()) {
    parser.value = std::make_unique<Kinetic<EngineType>>(system, event_filters);
  }
}

}  // namespace kinetic
}  // namespace clexmonte
}  // namespace CASM

#endif
