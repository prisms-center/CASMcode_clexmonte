# Changelog

All notable changes to `libcasm-clexmonte` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
