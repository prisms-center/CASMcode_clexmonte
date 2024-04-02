from collections.abc import Iterable
from typing import Optional, Union

import numpy as np

import libcasm.casmglobal as casmglobal
import libcasm.clexmonte as clexmonte
import libcasm.clexmonte.auto_configuration as autoconfig
import libcasm.configuration as casmconfig
from libcasm.monte import (
    RandomNumberEngine,
    ValueMap,
)
from libcasm.monte.events import (
    OccLocation,
)

from ._clexmonte_semigrand_canonical import (
    SemiGrandCanonicalCalculator as _SemiGrandCanonicalCalculator,
)
from ._clexmonte_semigrand_canonical import (
    SemiGrandCanonicalConditions,
    SemiGrandCanonicalPotential,
)


def make_initial_state(
    system: clexmonte.System,
    conditions: Union[dict, ValueMap, SemiGrandCanonicalConditions],
    dirs: str = "abc",
    min_volume: int = 1000,
    transformation_matrix_to_super: Optional[np.ndarray] = None,
    motif: Optional[dict] = None,
    configurations: Union[
        casmconfig.ConfigurationSet, Iterable[casmconfig.Configuration], None
    ] = None,
    tol: float = casmglobal.TOL,
):
    """Determine an appropriate initial state for Monte Carlo calculations

    Parameters
    ----------
    system: libcasm.clexmonte.System
        System data.
    conditions: Union[dict, ValueMap, SemiGrandCanonicalConditions]
        Semi-grand canonical ensemble thermodynamic conditions as a Python dict,
        :class:`ValueMap`, or :class:`SemiGrandCanonicalConditions`.
    dirs: str = "abc"
        The directions along which the initial supercell can be expanded ("a"
        corresponds to first supercell lattice vector, "b" the seoncd, and "c" the
        third).
    min_volume: int = 1000,
        The minimum volume of the results supercell, as integer multiples of the
        prim unit cell volume.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        If provided, only configurations that tile the corresponding supercell will
        be considered.
    motif: Optional[libcasm.configuration.Configuration] = None
        If not None, use the provided motif configuration as the initial state rather
        than finding the minimum potential configuration.
    configurations: Union[libcasm.configuration.ConfigurationSet, \
    Iterable[libcasm.configuration.Configuration], None] = None
        If not None, the candidate motif configurations which are searched for the
        minimum potential configuration. Must be a
        :class:`~libcasm.configuration.ConfigurationSet` or an iterable of
        :class:`~libcasm.configuration.Configuration`. If multiple approximately
        equal minimum potential configurations are found, the first encountered is
        used. If None, the default configuration is used.
    tol: float :data:`~libcasm.casmglobal.TOL`
        Tolerance for identifying configurations with approximately equal potential.

    Returns
    -------
    initial_state: clexmonte.MonteCarloState
        Initial Monte Carlo state according to the specified parameters.
    """
    potential = SemiGrandCanonicalPotential(system)
    if isinstance(conditions, dict):
        conditions = ValueMap(conditions)
    if isinstance(conditions, ValueMap):
        sgc_conditions = SemiGrandCanonicalConditions(
            composition_converter=system.composition_converter,
        )
        sgc_conditions.set_all(conditions, is_increment=False)
        conditions = sgc_conditions
    if not isinstance(conditions, SemiGrandCanonicalConditions):
        raise Exception("Invalid conditions type")

    return autoconfig.make_initial_state(
        system=system,
        potential=potential,
        conditions=conditions,
        dirs=dirs,
        min_volume=min_volume,
        transformation_matrix_to_super=transformation_matrix_to_super,
        motif=motif,
        configurations=configurations,
        tol=tol,
    )


class SemiGrandCanonicalCalculator(_SemiGrandCanonicalCalculator):
    def run(
        self,
        state: clexmonte.MonteCarloState,
        run_manager: clexmonte.RunManager,
        occ_location: Optional[OccLocation] = None,
    ):
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
        state: clexmonte.MonteCarloState,
        sampling_fixture_params: clexmonte.SamplingFixtureParams,
        engine: Optional[RandomNumberEngine] = None,
        occ_location: Optional[OccLocation] = None,
    ):
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
        run_manager = clexmonte.RunManager(
            engine=engine,
            sampling_fixture_params=[sampling_fixture_params],
        )
        run_manager = self.run(
            state=state,
            run_manager=run_manager,
            occ_location=occ_location,
        )
        return run_manager.sampling_fixtures[0]
