Overview
========

The libcasm-clexmonte package provides the CASM cluster expansion Monte Carlo
implementations.

General approach
----------------

1. Construct a :class:`~libcasm.clexmonte.System`.

   The System class holds data defining the crystal system and its properties, independent of a particular Monte Carlo method or particular state of the sytem. This includes:

   - the prim (:class:`~libcasm.configuration.Prim`),
   - the composition axes (:class:`~libcasm.composition.CompositionConverter`),
   - basis sets (:class:`~libcasm.clexulator.Clexulator`),
   - correlations calculators (:class:`~libcasm.clexulator.Correlations`),
   - cluster expansion calculators (:class:`~libcasm.clexulator.ClusterExpansion`), and
   - order parameter calculators (:class:`~libcasm.clexulator.OrderParameter`).

   The system data also includes other data such as local basis sets, local cluster expansions, and kinetic Monte Carlo events.

2. Construct a :class:`~libcasm.clexmonte.MonteCalculator`.

   The MonteCalculator class provides a standardized interface to Monte Carlo method
   implementations and all the data needed for the implementation. This includes:

   - the system data (:class:`~libcasm.clexmonte.System`),
   - the state data for the current state
     (:class:`~libcasm.clexmonte.StateData`),
   - the potential calculator (:class:`~libcasm.clexmonte.MontePotential`),
   - the standard sampling functions
     (:class:`~libcasm.monte.sampling.StateSamplingFunction` and
     :class:`~libcasm.monte.sampling.jsonStateSamplingFunction`) provided by the
     implementation,
   - results analysis functions (:class:`~libcasm.clexmonte.ResultsAnalysisFunction`), and
   - functions to run the Monte Carlo method (:func:`~libcasm.clexmonte.MonteCalculator.run_fixture` and :func:`~libcasm.clexmonte.MonteCalculator.run`).

3. Construct an initial :class:`~libcasm.clexmonte.MonteCarloState`.

   The MonteCarloState data structure combines:

   - a configuration (:class:`~libcasm.configuration.Configuration`), and
   - thermodynamic conditions (:class:`~libcasm.monte.ValueMap`).

   A MonteCarloState can be constructed by:

   - explicitly giving the exact configuration and conditions (using the :class:`~libcasm.clexmonte.MonteCarloState` constructor), or
   - using the :func:`~libcasm.clexmonte.make_initial_state` method to perform standard
     operations like finding the configuration with minimum potential, or fill a
     supercell with a certain shape or minimum volume with a motif configuration.

   For canonical and kinetic Monte Carlo calculations, it may be useful to:

   - use the :func:`~libcasm.clexmonte.enforce_composition` method to perturb an
     MonteCarloState configuration to match a desired composition, or
   - set the conditions of the MonteCarloState to match its configuration.

4. Construct one or more :class:`~libcasm.clexmonte.SamplingFixtureParams`

   A sampling fixture (:class:`~libcasm.clexmonte.SamplingFixture`) is used to sample
   data, store results, and check for completion during a Monte Carlo simulation.
   SamplingFixtureParams is a data structure that specifies all the parameters that
   control a :class:`~libcasm.clexmonte.SamplingFixture`. This includes:

   - sampling functions (:class:`~libcasm.monte.sampling.StateSamplingFunction` and
     :class:`~libcasm.monte.sampling.jsonStateSamplingFunction`), including both
     standard sampling functions provided by the implementation and user-provided
     custom sampling functions, which return the quantities (energy, composition, order
     parameters, etc.) sampled by the fixture,
   - sampling parameters (:class:`~libcasm.monte.sampling.SamplingParams`), specifying
     which sampling functions to evaluate and when the samples should be taken,
   - completion check parameters
     (:class:`~libcasm.monte.sampling.CompletionCheckParams`),
     which includes which sampled quantities should be converged, the requested
     absolute or relative precision level, how often to check, and minimum and maximum
     numbers of samples, steps or passes, computational or simulated time to run for,
   - results output, including where to write output files, whether to only write a
     summary with mean values and estimated precisions, or to also write all observations,
     or the trajectory of configurations at each sample time, and
   - status logging parameters, including whether, where, and how often to write a status
     log file with the most recent completion check results.

   In some cases it may be useful to use multiple sampling fixtures for a single Monte
   Carlo simulation. For instance, a sampling fixture for thermodynamic properties can
   be re-used and combined with a sampling fixture for kinetic properties during a
   kinetic Monte Carlo simulation, or sampling fixtures that sample different
   quantities at different intervals could be combined.

   The :class:`~libcasm.clexmonte.RunManager` class is used to holds one or more
   :class:`SamplingFixture` and make it easy to write a Monte Carlo method that does
   sampling and convergence checking according to each sampling fixture. A
   `global_cutoff` parameter determines if all sampling fixtures must complete for the
   Monte Carlo run to finish, or if the run should stop when any one sampling fixture
   completes.

   Additionally, the RunManager controls options for saving initial and final states of
   each run in order to enable re-starts of a series of runs and perform "dependent
   runs" where the final configuration of one run is used as the initial configuration
   of the next run at varying conditions.

5. Run one or more Monte Carlo simulations

   Simulations can be run at single state using:

   - :func:`~libcasm.clexmonte.MonteCalculator.run_fixture` when using a single sampling fixture, or
   - :func:`~libcasm.clexmonte.MonteCalculator.run` when using a RunManager.configuration

   Main results, the average value of sampled quantities and estimated precision,
   and the calculation of quantities like the heat capacity and susceptibility from
   fluctuations of energy and composition, are stored in a results summary file.
   The values calculated from each subsequent run are stored by appending to lists in a
   results summary file. Typically runs are organized along linear paths in
   thermodynamic conditions space (for instance increasing temperature at constant
   chemical potential), with one summary file for one linear path.

   This process can be automated by:

   - using the :func:`~libcasm.clexmonte.run_series` method to run a series of simulations along a path in conditions space,
   - using the `casm-flow <TODO>` package, which helps automate the process of setting
     up input files, submitting jobs to a cluster, and collecting, analyzing, and
     plotting results.

