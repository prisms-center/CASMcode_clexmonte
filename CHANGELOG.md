# Changelog

All notable changes to `libcasm-clexmonte` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.0] - 2025-08-14

### Changed

- Set pybind11~=3.0


## [2.1.0] - 2025-08-08

### Changed

- Build Linux wheels using manylinux_2_28 (previously manylinux2014)
- Removed Cirrus CI testing


## [2.0.0] - 2025-05-04

### Changed

- Build for Python 3.13
- Restrict requires-python to ">=3.9,<3.14"
- Run CI tests using Python 3.13
- Build MacOS arm64 wheels using MacOS 15
- Build Linux wheels using Ubuntu 24.04


## [2.0a6] - 2024-02-12

### Changed

- Removed erroneous import into `libcasm.clexmonte.misc.contexts`.


## [2.0a5] - 2024-02-11

### Added

- Added `verbose` argument to `System.from_dict` to print messages during parsing.
- Added `System.basis_set_name`, `System.local_basis_set_name`, `System.prototype_event`, and `libcasm.clexmonte.make_system_event_info` to get event info and support LocalConfiguration construction from SelectedEvent info.
- Added `SelectedEvent.to_dict` and `SelectedEvent.__repr__` methods.
- Added `MonteCalculator.event_info` attribute to store shared event info used for constructing and transforming LocalConfiguration.
- Added `MonteCalculator.make_local_configuration` to construct a LocalConfiguration from the current calculator state and current SelectedEvent.
- Added `EventStateCalculator` and `MonteEventData.set_custom_event_state_calculation`
- Added `"print_event_data_summary"` MonteCalculator parameter to select whether to print the event data summary.
- Added `libcasm.clexmonte.site_iterators` to use help use OccLocation more easily.
- Added options for handling abnormal events (events with no energy barrier between the initial and final states). This includes the `"abnormal_event_handling"` option and sub-options for the MonteCalculator constructor `params` argument, custom abnormal event handling functions for both encountered (calculated events that may not be selected) and selected events, writing and reading local configurations and event data for abnormal events (`MonteCalculator.read_abnormal_events`).
- Added a random number engine to MonteCalculator, to enable control of the random number generator used by the `"enforce.composition"` method in `make_canonical_initial_state`, and to use as the default for `MonteCalculator.run_fixture` if no engine is given explicitly. The MonteCalculator engine is replaced by the `run_manager` engine when `MonteCalculator.run` is called. 
- Added `libcasm.clexmonte.misc.contexts` to include a working directory context manager.
 
### Changed

- Changed MonteCalculator to a Python class that inherits from a C++ class, allowing for more flexibility in the Python interface.
- Moved `libcasm.clexmonte.print_selected_event_functions` to `MonteCalculator.print_selected_event_functions`.
- Changed `PrimEventData.sites` type to `SiteVector`
- Changed `PrimEventData.occ_init` and `PrimEventData.occ_final` type to `IntVector`
- Changed `EventState` attributes from `readonly` to `readwrite`.
- Changed `"kinetic"` MonteCalculator to not print the event data summary by default
- Enable surpressing "kinetic" MonteCalculator output by reducing verbosity to "quiet".

### Fixed

- Fix for unuseful error messages when `System.from_dict` failed to load a local basis set.
- Fixed setting `begin` parameter in `SamplingFixtureParams.sample_by_X` functions.


## [2.0a4] - 2024-12-17

### Fixed

- Fixed a bug that gave the wrong change in potential energy during semi-grand canonical Monte Carlo simulations. The bug was introduced in v2.0a3 only and is fixed in v2.0a4.
- Fixed parsing of other files listed in an input file, such as a coefficients file listed in a System input file. Now, errors and warnings for the file being parsed are properly shown.


## [2.0a3] - 2024-12-12

### Added

- Added "kinetic" MonteCalculator method for kinetic Monte Carlo simulations.
- Added setters to SamplingFixtureParams for setting sub-object parameters to make it easier to set parameters in a single line.
- Added selected event data sampling for KMC simulations.
- Added CASM::clexmonte::AllowedEventList and related classes so that all possible events do not need to be added to the KMC event selector. 
- Added CASM::clexmonte::EventDataSummary and libcasm.clexmonte.MonteEventDataSummary to summarize event data for KMC simulations.
- Optional "neighborlist" or "relative" impact table types.

### Changed

- The AllowedEventList method is made the default for the "kinetic" MonteCalculator method. The event data type can be selected using the `params` argument to the MonteCalculator constructor.
- Changed the enforce_composition method to avoid unnecessarily re-calcuating the composition at each step.


## [2.0a2] - 2024-07-17

### Fixed

- Updated for compatibility with libcasm-configuration 2.0a5



## [2.0a1] - 2024-07-17

This release creates the libcasm-clexmonte cluster expansion based Monte Carlo module. It includes:

- Canonical, semi-grand canonical, and kinetic Monte Carlo calculators
- Support for customizing potentials, including linear, quadratic, and correlation-matching terms 
- Metropolis and N-fold way implementations
- Support for customizing sampling and analysis functions

The distribution package libcasm-clexmonte contains a Python package (libcasm.clexmonte) that provides an interface to Monte Carlo simulation methods implemented in C++. The libcasm.clexmonte.MonteCalculator class currently provides access to simulations in the canonical and semi-grand canonical ensemble and will be expanded in the next releases to include additional methods.

This package may be installed via pip install, using scikit-build, CMake, and pybind11. This release also includes usage examples and API documentation, built using Sphinx.
