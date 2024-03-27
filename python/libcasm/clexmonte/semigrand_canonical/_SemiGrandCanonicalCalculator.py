from collections.abc import Iterable
from typing import Optional, Union

import numpy as np

import libcasm.casmglobal as casmglobal
import libcasm.clexmonte as clexmonte
import libcasm.clexmonte.auto_configuration as autoconfig
import libcasm.configuration as casmconfig
from libcasm.monte import (
    ValueMap,
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
        conditions = ValueMap.from_dict(conditions)
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
    def make_initial_state(
        self,
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
            If not None, use the provided motif configuration as the initial state
            rather than finding the minimum potential configuration.
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
        return make_initial_state(
            system=self.system,
            conditions=conditions,
            dirs=dirs,
            min_volume=min_volume,
            transformation_matrix_to_super=transformation_matrix_to_super,
            motif=motif,
            configurations=configurations,
            tol=tol,
        )
