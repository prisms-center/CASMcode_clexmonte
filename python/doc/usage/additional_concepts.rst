Additional Concepts
===================

.. toctree::
    :maxdepth: 2
    :hidden:

This page provides an overview of some additional concepts that may be useful for those creating custom Monte Carlo methods with CASM.


Customized calculation types
----------------------------

A custom calculation type, for instance implementing a new potential, can be constructed by copying an existing calculation type directory and modifying the potential, event definitions, etc. If it is of general interest, a pull request may be accpted to merge the new calculation type back into CASM_clexmonte.


OccLocation
-----------
Internally, the `OccLocation` class is used to track occupants in a configuration. This allows for efficient selection of events in canonical Monte Carlo calculations and for tracking displacement of atoms including periodic boundary crossings in kinetic Monte Carlo calculations. An `OccLocation` instance should be constructed and initialized based on the initial configuration before each Monte Carlo run. Configuration changes for all accepted swaps or events are then applied through `OccLocation::apply` so that internal data structures can be updated to reflect the changes.


Run manager
-----------

A `RunManager` instance holds all sample fixtures and applies them during a calculation. It continues a run until all sampling fixtures are complete. It also stores a list of final states.


State generators
----------------

A state generator is used to construct the initial state for each run in a series from user provided parameters and results from previous calculation runs. The simplest state generator, `IncrementalConditionsStateGenerator`, takes parameters specifying an initial configuration, including supercell size and shape and degree of freedom (DoF) values, initial thermodynamic conditions, incremental conditions, and the number of total states to be generated. Each generated state is used as the input to a new Monte Carlo calculation run.

In the future, additional state generators will be added.


Customized conditions
---------------------

Currently a single `Conditions` class is used for all types of ensembles, ignoring parameters that are not applicable. The `Conditions` class may be converted to / from a `ValueMap` data structure which is a mapping of name to boolean, scalar, vector, or matrix values. For compatibility with `IncrementalConditionsStateGenerator` and results output, when a new condition type is needed the `Conditions` class and conversion methods should be updated.


Variable conditions
-------------------

A future update to allow variable conditions is planned. It will...
