KMC customization
=================

Sampling functions
------------------

Sampling functions take individual observations of the current configuration at each sampling time. Values of any dimension may be sampled by defining an unrolling convention and returning a vector of values. Each component is given a name for identifiation in the results output file. CASM can check for equilibration and convergence of each individual component. If a combined metric is of interest for convergence, a seperate sampling function should be constructed.

Examples of sampling functions...

Kinetic sampling functions have access to the current configuration and a limited amount of shared information about the state of the configuration at the last sample time: previous sample time, atom positions, etc. Additional data needed by a sampling function about the previous state may be saved in a std::shared_ptr.

Example of saving state...

To add a new sampling function...


Analysis functions
------------------

Analysis functions perform calculations using all collected observations once a Monte Carlo run is complete. Examples include calculation of measures of fluctuations such as the heat capacity and susceptibility.

Examples of analysis functions...

To add a new analysis function...
