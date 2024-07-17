# Changelog

All notable changes to `libcasm-clexmonte` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [2.0a1] - 2024-07-17

This release creates the libcasm-clexmonte cluster expansion based Monte Carlo module. It includes:

- Canonical, semi-grand canonical, and kinetic Monte Carlo calculators
- Support for customizing potentials, including linear, quadratic, and correlation-matching terms 
- Metropolis and N-fold way implementations
- Support for customizing sampling and analysis functions

The distribution package libcasm-clexmonte contains a Python package (libcasm.clexmonte) that provides an interface to Monte Carlo simulation methods implemented in C++. The libcasm.clexmonte.MonteCalculator class currently provides access to simulations in the canonical and semi-grand canonical ensemble and will be expanded in the next releases to include additional methods.

This package may be installed via pip install, using scikit-build, CMake, and pybind11. This release also includes usage examples and API documentation, built using Sphinx.
