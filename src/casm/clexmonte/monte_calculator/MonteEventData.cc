#include "casm/clexmonte/monte_calculator/MonteEventData.hh"

// -- Memory usage --
// https://stackoverflow.com/a/372525

#ifdef __linux__
#include <sys/sysinfo.h>
#endif

#ifdef __APPLE__
#include <mach/mach_init.h>
#include <mach/task.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#else
#include <sys/resource.h>
#endif

namespace CASM {

/// \brief The amount of memory currently being used by this process, in bytes.
///
/// By default, returns the full virtual arena, but if resident=true,
/// it will report just the resident set in RAM (if supported on that OS).
size_t memory_used(bool resident) {
#if defined(__linux__)
  // Ugh, getrusage doesn't work well on Linux.  Try grabbing info
  // directly from the /proc pseudo-filesystem.  Reading from
  // /proc/self/statm gives info on your own process, as one line of
  // numbers that are: virtual mem program size, resident set size,
  // shared pages, text/code, data/stack, library, dirty pages.  The
  // mem sizes should all be multiplied by the page size.
  size_t size = 0;
  FILE *file = fopen("/proc/self/statm", "r");
  if (file) {
    unsigned long vm = 0;
    fscanf(file, "%ul", &vm);  // Just need the first num: vm size
    fclose(file);
    size = (size_t)vm * getpagesize();
  }
  return size;

#elif defined(__APPLE__)
  // Inspired by:
  // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info,
            &t_info_count);
  size_t size = (resident ? t_info.resident_size : t_info.virtual_size);
  return size;

#elif defined(_WINDOWS)
  // According to MSDN...
  PROCESS_MEMORY_COUNTERS counters;
  if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)))
    return counters.PagefileUsage;
  else
    return 0;

#else
  // No idea what platform this is
  return 0;  // Punt
#endif
}

std::string convert_size(size_t size_bytes) {
  // Based on # https://stackoverflow.com/a/14822210
  if (size_bytes == 0) {
    return "0B";
  }
  std::vector<std::string> size_name = {"B",   "KiB", "MiB", "GiB", "TiB",
                                        "PiB", "EiB", "ZiB", "YiB"};
  double i = std::floor(std::log2(size_bytes) / std::log2(1024));
  double s = size_bytes / std::pow(1024, i);
  std::stringstream ss;
  ss << std::fixed << std::setprecision(2) << s
     << size_name[static_cast<int>(i)];
  return ss.str();
}

}  // namespace CASM

// -- End Memory usage --

namespace CASM {
namespace clexmonte {

EventTypeStats::EventTypeStats(
    std::vector<std::string> const &_partion_names_by_type,
    std::vector<std::string> const &_partion_names_by_equivalent_index,
    double _initial_begin, double _bin_width, bool _is_log, Index _max_size)
    : n_total(0),
      min(0.0),
      max(0.0),
      sum(0.0),
      mean(0.0),
      hist_by_type(_partion_names_by_type, _initial_begin, _bin_width, _is_log,
                   _max_size),
      hist_by_equivalent_index(_partion_names_by_equivalent_index,
                               _initial_begin, _bin_width, _is_log, _max_size) {
}

void EventTypeStats::insert(int partition_by_type,
                            int partition_by_equivalent_index, double value) {
  n_total += 1;
  min = std::min(min, value);
  max = std::max(max, value);
  sum += value;
  mean = sum / n_total;
  hist_by_type.insert(partition_by_type, value);
  hist_by_equivalent_index.insert(partition_by_equivalent_index, value);
}

EventDataSummary::EventDataSummary(
    std::shared_ptr<StateData> const &_state_data,
    MonteEventData const &_event_data, double energy_bin_width,
    double freq_bin_width, double rate_bin_width)
    : state_data(_state_data),
      event_data(_event_data),
      prim_event_list(event_data.prim_event_list()) {
  for (Index i = 0; i < prim_event_list.size(); ++i) {
    auto type = type_key(i);
    auto equiv = equiv_key(i);

    all_types.insert(type);
    equiv_keys_by_type[type].insert(equiv);
    all_equiv_keys.insert(equiv);

    n_possible.by_type[type] = 0.0;
    n_possible.by_equivalent_index[equiv] = 0.0;

    n_allowed.by_type[type] = 0.0;
    n_allowed.by_equivalent_index[equiv] = 0.0;

    n_not_normal.by_type[type] = 0.0;
    n_not_normal.by_equivalent_index[equiv] = 0.0;

    rate.by_type[type] = 0.0;
    rate.by_equivalent_index[equiv] = 0.0;

    n_impact.by_type[type] = 0.0;
    n_impact.by_equivalent_index[equiv] = 0.0;
  }

  n_events_allowed = 0.0;
  n_events_possible = 0.0;
  event_list_size = event_data.event_list().size();
  total_rate = event_data.event_list().total_rate();
  mean_time_increment = 1.0 / total_rate;

  resident_bytes_used = memory_used(true);
  virtual_bytes_used = memory_used(false);

  SelectedEventInfo info(event_data.prim_event_list());
  info.make_indices_by_type();
  to_event_type = *info.prim_event_index_to_index;
  event_type_names = info.partition_names;

  info.make_indices_by_equivalent_index();
  to_equivalent_index = *info.prim_event_index_to_index;
  equivalent_index_names = info.partition_names;

  stats_labels.push_back("dE_final");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("dE_activated");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("Ekra");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, energy_bin_width /* bin_width */,
                     false /* is_log */);

  stats_labels.push_back("freq");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, freq_bin_width /* bin_width */,
                     true /* is_log */);

  stats_labels.push_back("rate");
  stats.emplace_back(event_type_names, equivalent_index_names,
                     0.0 /* initial_begin */, rate_bin_width /* bin_width */,
                     true /* is_log */);

  for (EventID const &id : event_data.event_list()) {
    EventState const &state = event_data.event_state(id);
    _add_count(id, state);
    _add_impact(id, state);
    _add_stats(id, state);
  }
}

void EventDataSummary::_add_count(EventID const &id, EventState const &state) {
  auto type = type_key(id);
  auto equiv = equiv_key(id);

  // double increment = 1.0 / state_data->n_unitcells;
  Index increment = 1;

  if (prim_event_list[id.prim_event_index].is_forward) {
    n_possible.by_type[type] += increment;
    n_possible.by_equivalent_index[equiv] += increment;
    n_events_possible += increment;
  }
  if (state.is_allowed) {
    n_events_allowed += increment;
    n_allowed.by_type[type] += increment;
    n_allowed.by_equivalent_index[equiv] += increment;
    if (!state.is_normal) {
      n_not_normal.by_type[type] += increment;
      n_not_normal.by_equivalent_index[equiv] += increment;
    }
  }
  rate.by_type[type] += state.rate;
  rate.by_equivalent_index[equiv] += state.rate;
}

void EventDataSummary::_add_impact(EventID const &id, EventState const &state) {
  TypeKey type = type_key(id);
  EquivKey equiv = equiv_key(id);

  // double increment = 1.0 / state_data->n_unitcells;
  Index increment = 1;

  // -- impact_table.by_type --
  for (EventID const &impact_id : event_data.event_impact(id)) {
    TypeKey impact_type = type_key(impact_id);
    n_impact.by_type[type] += increment;
    impact_table.by_type[type][impact_type] += increment;
  }

  // -- impact_table.by_equivalent_index --
  for (EventID const &impact_id : event_data.event_impact(id)) {
    EquivKey impact_equiv = equiv_key(impact_id);
    n_impact.by_equivalent_index[equiv] += increment;
    impact_table.by_equivalent_index[equiv][impact_equiv] += increment;
  }

  // -- neighborhood_size_total --
  if (neighborhood_size_total.count(type) == 0) {
    PrimEventData const &prim_event_data =
        event_data.prim_event_list()[id.prim_event_index];
    System const &system = *event_data.system();
    OccEventTypeData const &event_type_data =
        get_event_type_data(system, prim_event_data.event_type_name);
    clust::IntegralCluster phenom = make_cluster(prim_event_data.event);

    std::set<xtal::UnitCellCoord> total_hood;
    ClexData const &clex_data = get_clex_data(system, "formation_energy");
    expand(phenom, total_hood, *clex_data.cluster_info, clex_data.coefficients);
    neighborhood_size_formation_energy[type] = total_hood.size();

    // local basis set dependence
    LocalMultiClexData const &local_multiclex_data =
        get_local_multiclex_data(system, event_type_data.local_multiclex_name);
    std::set<xtal::UnitCellCoord> kra_hood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index, "kra");
    neighborhood_size_kra[type] = kra_hood.size();
    total_hood.insert(kra_hood.begin(), kra_hood.end());

    std::set<xtal::UnitCellCoord> freq_hood = get_required_update_neighborhood(
        system, local_multiclex_data, prim_event_data.equivalent_index, "freq");
    neighborhood_size_freq[type] = freq_hood.size();
    total_hood.insert(freq_hood.begin(), freq_hood.end());

    neighborhood_size_total[type] = total_hood.size();
  }
}

void EventDataSummary::_add_stats(EventID const &id, EventState const &state) {
  if (!state.is_allowed) {
    return;
  }

  int t = to_event_type[id.prim_event_index];
  int e = to_equivalent_index[id.prim_event_index];

  // order determined by constructor
  int i = 0;
  stats[i++].insert(t, e, state.dE_final);
  stats[i++].insert(t, e, state.dE_activated);
  stats[i++].insert(t, e, state.Ekra);
  stats[i++].insert(t, e, state.freq);
  stats[i++].insert(t, e, state.rate);
}

void print(std::ostream &out, EventDataSummary const &event_data_summary) {
  EventDataSummary const &x = event_data_summary;
  out << "Event data summary:\n";
  out << "- Number of unitcells = " << x.state_data->n_unitcells << std::endl;
  // Number of events
  out << "- Number of events (total) = " << x.n_events_allowed << std::endl;
  out << "- Number of events (by type):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = "
        << x.n_allowed.by_type.at(event_type_name);
    if (x.n_not_normal.by_type.at(event_type_name) != 0.0) {
      out << " (** not_normal = " << x.n_not_normal.by_type.at(event_type_name)
          << " **)";
    }
    out << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      out << "    - " << event_type_name << "." << equivalent_index << " = "
          << x.n_allowed.by_equivalent_index.at(equiv_key);
      if (x.n_not_normal.by_equivalent_index.at(equiv_key) != 0.0) {
        out << " (** not_normal = "
            << x.n_not_normal.by_equivalent_index.at(equiv_key) << " **)";
      }
      out << std::endl;
    }
  }
  // Rate info:
  out << "- Event rate (total) (1/s) = " << x.total_rate << std::endl;
  out << "- Event rate (by type) (1/s):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = "
        << x.rate.by_type.at(event_type_name) << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      out << "    - " << event_type_name << "." << equivalent_index << " = "
          << x.rate.by_equivalent_index.at(equiv_key) << std::endl;
    }
  }
  out << "- Mean time increment (total) (s) = " << x.mean_time_increment
      << std::endl;
  out << "- Mean time increment (by type) (s):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = "
        << 1.0 / x.rate.by_type.at(event_type_name) << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      out << "    - " << event_type_name << "." << equivalent_index << " = "
          << 1.0 / x.rate.by_equivalent_index.at(equiv_key) << std::endl;
    }
  }
  // Memory usage:
  out << "- Memory used (RAM) = " << convert_size(x.resident_bytes_used)
      << std::endl;
  out << "- Event list size = " << x.event_list_size << std::endl;
  // Impact neighborhood sizes:
  out << "- Impact neighborhood sizes (#sites): total (Ef / Ekra / "
         "freq)"
      << std::endl;
  for (auto const &pair : x.neighborhood_size_total) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = " << pair.second << " ("
        << x.neighborhood_size_formation_energy.at(event_type_name) << " / "
        << x.neighborhood_size_kra.at(event_type_name) << " / "
        << x.neighborhood_size_freq.at(event_type_name) << ")" << std::endl;
  }
  // Impact list:
  out << "- Impact list sizes (total):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = "
        << x.n_impact.by_type.at(event_type_name) << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      out << "    - " << event_type_name << "." << equivalent_index << " = "
          << x.n_impact.by_equivalent_index.at(equiv_key) << std::endl;
    }
  }
  // Impact table:
  out << "- Impact list sizes (by type):" << std::endl;
  for (auto const &pair : x.equiv_keys_by_type) {
    std::string event_type_name = pair.first;
    out << "  - " << event_type_name << " = ";
    for (auto const &pair2 : x.impact_table.by_type.at(event_type_name)) {
      out << pair2.second << " ";
    }
    out << std::endl;

    for (auto const &equiv_key : pair.second) {
      Index equivalent_index = equiv_key.second;
      out << "    - " << event_type_name << "." << equivalent_index << " = ";
      for (auto const &pair2 :
           x.impact_table.by_equivalent_index.at(equiv_key)) {
        out << pair2.second << " ";
      }
      out << std::endl;
    }
  }
  out << std::endl;
}

}  // namespace clexmonte
}  // namespace CASM