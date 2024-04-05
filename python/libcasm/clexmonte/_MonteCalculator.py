from typing import Optional

from libcasm.monte import (
    RandomNumberEngine,
)
from libcasm.monte.events import (
    OccLocation,
)

from ._clexmonte_monte_calculator import (
    MonteCalculator as _MonteCalculator,
)
from ._clexmonte_run_management import (
    RunManager,
    SamplingFixture,
    SamplingFixtureParams,
)
from ._clexmonte_state import (
    MonteCarloState,
)


class MonteCalculator(_MonteCalculator):
    """Interface for running Monte Carlo calculations

    The MonteCalculator class provides a common interface for different
    Monte Carlo method implementations.
    """

    def run(
        self,
        state: MonteCarloState,
        run_manager: RunManager,
        occ_location: Optional[OccLocation] = None,
    ) -> RunManager:
        """Perform a single run, evolving the input state

        Parameters
        ----------
        state : libcasm.clexmonte.MonteCarloState
            The input state.
        run_manager: libcasm.clexmonte.RunManager
            Specifies sampling and convergence criteria and collects results
        occ_location: Optional[libcasm.monte.events.OccLocation] = None
            Current occupant location list. If provided, the user is
            responsible for ensuring it is up-to-date with the current
            occupation of `state`. It is used and updated during the run.
            If None, an occupant location list is generated for the run.

        Returns
        -------
        run_manager: libcasm.clexmonte.RunManager
            The input `run_manager` with collected results.
        """
        return self._run(
            state=state, run_manager=run_manager, occ_location=occ_location
        )

    def run_fixture(
        self,
        state: MonteCarloState,
        sampling_fixture_params: SamplingFixtureParams,
        engine: Optional[RandomNumberEngine] = None,
        occ_location: Optional[OccLocation] = None,
    ) -> SamplingFixture:
        """Perform a single run, with a single sampling fixture, evolving the \
        input state

        Parameters
        ----------
        state : libcasm.clexmonte.MonteCarloState
            The input state.
        sampling_fixture_params: libcasm.clexmonte.SamplingFixtureParams
            Specifies sampling and convergence criteria and collects results.
        engine: Optional[libcasm.monte.RandomNumberEngine] = None
            Optional random number engine to use. If None, one is constructed and
            seeded from std::random_device.
        occ_location: Optional[libcasm.monte.events.OccLocation] = None
            Current occupant location list. If provided, the user is
            responsible for ensuring it is up-to-date with the current
            occupation of `state`. It is used and updated during the run.
            If None, an occupant location list is generated for the run.

        Returns
        -------
        sampling_fixture: libcasm.clexmonte.SamplingFixture
            A SamplingFixture with collected results.

        """
        if engine is None:
            engine = RandomNumberEngine()
        run_manager = RunManager(
            engine=engine,
            sampling_fixture_params=[sampling_fixture_params],
        )
        run_manager = self.run(
            state=state,
            run_manager=run_manager,
            occ_location=occ_location,
        )
        return run_manager.sampling_fixtures[0]
