# MonteCalculatorType = NewType("MonteCalculator", None)


# class StateData(_StateData):
#     """Access state-specific data used in a Monte Carlo method"""
#
#     def __init__(
#         self,
#         system: System,
#         state: MonteCarloState,
#         occ_location: Optional[OccLocation] = None,
#         update_atoms: bool = False,
#     ):
#         """
#
#         .. rubric:: Constructor
#
#         Parameters
#         ----------
#         system : libcasm.clexmonte.System
#             Cluster expansion model system data.
#         state : libcasm.clexmonte.MonteCarloState
#             The input state.
#         occ_location: Optional[libcasm.monte.events.OccLocation] = None
#               Current occupant location list. If provided, the user is
#               responsible for ensuring it is up-to-date with the current
#               occupation of `state` and it is used and updated during the run.
#               If None, no occupant location list is stored. The occupant
#               location list is not required for evaluating the potential.
#         update_atoms : bool = False
#             Use `True` in kinetic Monte Carlo calculations to track species
#             location changes. Otherwise use `False`.
#         """
#         super().__init__(
#             system=system,
#             state=state,
#             occ_location=occ_location,
#             update_atoms=update_atoms,
#         )


# class MontePotential(_MontePotential):
#     """Interface to potential calculators
#
#     The MontePotential class provides a common interface for different
#     Monte Carlo potential calculating implementations.
#     """
#
#     def __init__(
#         self,
#         calculator: MonteCalculatorType,
#         state: MonteCarloState,
#     ):
#         """
#
#         .. rubric:: Constructor
#
#         Parameters
#         ----------
#         calculator : libcasm.clexmonte.MonteCalculator
#             Monte Carlo calculator which implements the potential.
#         state : libcasm.clexmonte.MonteCarloState
#             The state to be calculated.
#         """
#         super().__init__(calculator=calculator, state=state)

#
# class MonteCalculator(_MonteCalculator):
#     """Interface for running Monte Carlo calculations
#
#     The MonteCalculator class provides a common interface for different
#     Monte Carlo method implementations.
#     """
#
#     # def __init__(
#     #     self,
#     #     method: str,
#     #     system: System,
#     #     params: Optional[dict] = None,
#     # ):
#     #     """
#     #     .. rubric:: Constructor
#     #
#     #     Parameters
#     #     ----------
#     #     method : str
#     #         Monte Carlo method name. The options are:
#     #
#     #         - "semigrand_canonical": `Semi-grand canonical ensemble <todo>`_
#     #         - "canonical": `Canonical ensemble <todo>`_
#     #         - TODO "kinetic": `Kinetic Monte Carlo <todo>`_
#     #         - TODO "flex": Allows a range of custom potentials, including
#     #           composition and order parameter variance-constrained potentials,
#     #           and correlation-matching potentials
#     #
#     #     system : libcasm.clexmonte.System
#     #         Cluster expansion model system data. The required data depends on
#     #         the calculation method. See links under `method` for what system
#     #         data is required for each method.
#     #
#     #     params: Optional[dict] = None
#     #         Monte Carlo calculation method parameters. Expected values
#     #         depends on the calculation method. Options, with links to
#     #         parameter documentation and examples, include:
#     #
#     #         - "enumeration": `Save states <todo>`_ encountered during the
#     #           calculation.
#     #     """
#     #     super().__init__(method=method, system=system, params=params)
#
#     # @property
#     # def state_data(self):
#     #     """StateData : The current state data."""
#     #     return StateData.self._state_data
#     #
#     # @property
#     # def potential(self):
#     #     """MontePotential : The potential calculator for the current state."""
#     #     return self._potential
#     #
#     # def run(
#     #     self,
#     #     state: MonteCarloState,
#     #     run_manager: RunManager,
#     #     occ_location: Optional[OccLocation] = None,
#     # ) -> RunManager:
#     #     """Perform a single run, evolving the input state
#     #
#     #     Parameters
#     #     ----------
#     #     state : libcasm.clexmonte.MonteCarloState
#     #         The input state.
#     #     run_manager: libcasm.clexmonte.RunManager
#     #         Specifies sampling and convergence criteria and collects results
#     #     occ_location: Optional[libcasm.monte.events.OccLocation] = None
#     #         Current occupant location list. If provided, the user is
#     #         responsible for ensuring it is up-to-date with the current
#     #         occupation of `state`. It is used and updated during the run.
#     #         If None, an occupant location list is generated for the run.
#     #
#     #     Returns
#     #     -------
#     #     run_manager: libcasm.clexmonte.RunManager
#     #         The input `run_manager` with collected results.
#     #     """
#     #     return self._run(
#     #         state=state, run_manager=run_manager, occ_location=occ_location
#     #     )
#     #
#     # def run_fixture(
#     #     self,
#     #     state: MonteCarloState,
#     #     sampling_fixture_params: SamplingFixtureParams,
#     #     engine: Optional[RandomNumberEngine] = None,
#     #     occ_location: Optional[OccLocation] = None,
#     # ) -> SamplingFixture:
#     #     """Perform a single run, with a single sampling fixture, evolving the \
#     #     input state
#     #
#     #     Parameters
#     #     ----------
#     #     state : libcasm.clexmonte.MonteCarloState
#     #         The input state.
#     #     sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
#     #         Specifies sampling and convergence criteria and collects results.
#     #     engine: Optional[libcasm.monte.RandomNumberEngine] = None
#     #         Optional random number engine to use. If None, one is constructed and
#     #         seeded from std::random_device.
#     #     occ_location: Optional[libcasm.monte.events.OccLocation] = None
#     #         Current occupant location list. If provided, the user is
#     #         responsible for ensuring it is up-to-date with the current
#     #         occupation of `state`. It is used and updated during the run.
#     #         If None, an occupant location list is generated for the run.
#     #
#     #     Returns
#     #     -------
#     #     sampling_fixture: libcasm.clexmonte.SamplingFixture
#     #         A SamplingFixture with collected results.
#     #
#     #     """
#     #     if engine is None:
#     #         engine = RandomNumberEngine()
#     #     run_manager = RunManager(
#     #         engine=engine,
#     #         sampling_fixture_params=[sampling_fixture_params],
#     #     )
#     #     run_manager = self.run(
#     #         state=state,
#     #         run_manager=run_manager,
#     #         occ_location=occ_location,
#     #     )
#     #     return run_manager.sampling_fixtures[0]
