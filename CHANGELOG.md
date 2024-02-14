# Changelog

All notable changes to `libcasm-clexmonte` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed

- Fixed CASM::config::make_distinct_cluster_sites

### Changed

- Changed libcasm.configuration.Prim.is_atomic, which was returning has_anisotropic_occupants
- Changed libcasm.configuration.is_canonical_configuration arguments for providing a subgroup

### Added

- Added prototype and include_subclusters attributes to libcasm.clusterography.IntegralClusterOrbitGenerator
- Added ConfigEnumAllOccupations, ConfigEnumInfo, ScelEnum, make_all_distinct_periodic_perturbations, and make_distinct_cluster_sites to libcasm.enumerate
- Added make_all_super_configurations_by_subsets and make_distinct_super_configurations to libcasm.configuration


## [v2.0a2] - 2023-12-11

### Fixed

- Fix bug in irrep decomposition
- Fix comparison of Configuration with equivalent but distinct supercell

### Added

- Added more irrep decomposition, DoF space analysis, and config space analysis tests
- Added options to config_space_analysis and dof_space_analysis methods to specify default occupation mode on a site or sublattice basis
- Added CASM::config::make_dof_space_rep and libcasm.configuration.make_dof_space_rep
- Added libcasm.configuration.ConfigurationWithProperties, and methods for libcasm.configuration.SupercellSymOp to act on ConfigurationWithProperties
- Added site_index_converter and unitcell_index_converter accessors to libcasm.Supercell.
- Added to_structure and from_structure to libcasm.configuration.Configuration and libcasm.configuration.ConfigurationWithProperties for conversions between atomic structures and configuration
- Added more access to matrix reps from PrimSymInfo in libcasm.configuration.Prim
- Added to_index_list, to_index_set, sort, sorted, is_sorted, __rmul__ to libcasm.clusterography.Cluster
- Added make_periodic_equivalence_map, make_periodic_equivalence_map_indices, make_local_equivalence_map, and make_local_equivalence_map_indices to libcasm.clusterography
- Added to_dict and from_dict methods to libcasm.configuration.Prim

### Changed

- Changed libcasm.clusterography.make_prim_periodic_orbits to make_periodic_orbits

### Deprecated

- Deprecated to_json and from_json methods of libcasm.configuration.Prim


## [v2.0a1] - Unreleased

This release creates the libcasm-clexmonte cluster expansion based Monte Carlo module. It includes:

- Canonical, semi-grand canonical, and kinetic Monte Carlo calculators
- Support for customizing potentials, including linear, quadratic, and correlation-matching terms 
- Metropolis and N-fold way implementations
- Support for customing sampling and analysis functions

The distribution package libcasm-clexmonte contains several Python packages of use for configuration comparison and enumeration:

- libcasm.canonical
- libcasm.semigrandcanonical
- libcasm.kinetic
- 
- libcasm.custom
- libcasm.enumerate

This package may be installed via pip install, using scikit-build, CMake, and pybind11. This release also includes usage examples and API documentation, built using Sphinx.
