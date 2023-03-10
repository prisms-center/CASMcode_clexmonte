Concepts
========

.. toctree::
    :maxdepth: 2
    :hidden:

    additional_concepts

This page provides an overview of concepts for CASM users.


Calculation Type
----------------

Examples:

- :class:`libcasm.canonical.Canonical`
- :class:`libcasm.semi_grand_canonical.SemiGrandCanonical`
- :class:`libcasm.kinetic.Kinetic`

A calculation type:

- Stores the data needed for a particular type of calculation
- Provides standard sampling, analysis, and modifying functions
- Implements a `run` method to run a calculation at a single set of thermodynamic conditions
- Implements a `run_series` method to multiple calculations at a series of thermodynamic conditions


Configuration
-------------

A configuration is specified by a prim, supercell, and the values of all degrees of freedom (DoF). This includes all local DoF at sites in that supercell and all global DoF.


Conditions
----------

Thermodynamic conditions that are being controlled in the current ensemble. For example:

- In the canonical ensemble, the conditions being controlled are temperature and composition
- In the semi-grand canonical ensemble, the conditions being controlled are temperature and exchange chemical potential
- In kinetic Monte Carlo calculations, the conditions being controlled are typically temperature and composition
- In a variance-constrained ensemble, the conditions also include the bias potential curvature and position


State
-----

A state is used to refer to a configuration and the applied thermodynamic conditions. Additional properties, such as the formation energy, may also be associated with a state.


Sampling Function
-----------------

Sampling functions take individual observations of the current state each time a sample is requested. Values of any dimension may be sampled by defining an unrolling convention and returning a vector of values. Each component is given a name for identifiation in the results output file. CASM can check for equilibration and convergence of each individual component. If a combined metric is of interest for convergence, a seperate sampling function should be constructed.


Analysis functions
------------------

Analysis functions perform calculations using all collected observations once a Monte Carlo run is complete. Examples include calculation of measures of fluctuations such as the heat capacity and susceptibility.


Configuration modifying functions
---------------------------------

Configuration modifying functions transform a configuration based on the current configuration and thermodynamic conditions. Currently this is limited to the function `set_mol_composition` which is applied before canonical and kinetic Monte Carlo calculations to modify the initial state such that the composition matches the thermodynamic conditions. Additional modifying functions may be added in the future.


Swaps
-----

For canonical and semi-grand canonical Monte Carlo calculations, CASM generates a list of allowed swaps. Currently swaps are limited to changing the occupation on one or two sites. A swap is defined by two "candidates" for swapping, which are specified using:

- Asymmetric unit index
- Species index

In typical usage, the species index is an index into a list of the types of atoms allowed by the prim. More generally, it is an index into a list of species that cannot be translated onto each other (for example, to distinguish particles with different spin or different orientations of a molecule).

For canonical swaps, which are used to swap species on two different sites, the candidate species must be different and the must be allowed on both asymmetric unit sites (if they differ).

For semi-grand canonical swaps, which are used to swap the species on a single site, the candidate species must be different and the asymmetric unit must be the same.


Events
------

For kinetic Monte Carlo and N-fold way calculations, CASM can generate a list of allowed events. Each event consists of a list of trajectories defining which occupants are moving, which sites they start on, and which sites they end on. Sites are defined using integral site coordinates, `(i,j,k,b)`, where `i,j,k` specify a unit cell as multiples of the primitive lattice vectors, and `b` is the sublattice index. A typical example is:

- A vacancy B-atom exchange event (Va-B) is defined as:
  - Va on site (0,0,0,0) moves to be Va on site (1,0,0,0)
  - B on site (1,0,0,0) moves to be B on site (0,0,0,0)

If a system allows A-Va and B-Va exchange events, then two OccEvent are used.

This definition more generally allows specifying the particular occupant on each site before and after the event and which atoms move to which position. This allows specifying more complicated events such as multi-occupant hops (i.e. triplet exchange), molecule re-orientation, and dumbbell breakup or formation. Moves to or from a resevoir may also be specified.

This specification does not specify the pathway, so if for instance there are two pathways connecting the same end states, then either two events should be included or the collective effect parameterized for a single event. Additionally, symmetry is determined by the end states only, so if the pathway breaks symmetry that must be accounted for specially.


Sampling fixtures
-----------------

A sampling fixture is collection of sampling functions, parameters controlling the frequency at which samples should be taken, calculation completion criteria, and analysis functions to be applied when a run is completed. Completion criteria may include minimum and maximum numbers of steps, passes, and samples or simulation time. They may also include convergence criteria in terms of a requested precision. Sampling fixtures also specify how output should be written and whether and how frequently a status file should be output during calculations.

Multiple sampling fixtures may be used in a calculation. This may be used for convenience to re-use common input files, or to specify different frequencies for sampling. Sampling fixtures which reach their own completion criteria continue to collect samples and converge results while other fixtures are not complete.
