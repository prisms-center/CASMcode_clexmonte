"""Find minimum potential configurations"""
from collections.abc import Iterable
from typing import Optional, Union

import numpy as np

import libcasm.casmglobal as casmglobal
import libcasm.configuration as casmconfig
import libcasm.monte as monte

from ._clexmonte_state import (
    MonteCarloState,
)
from ._MonteCalculator import MonteCalculator


def scale_supercell(
    transformation_matrix_to_super: np.ndarray,
    dirs: str = "abc",
    min_volume: int = 1000,
) -> np.ndarray:
    """Scale supercell along specified axes to have specified minimum volume

    Parameters
    ----------
    transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
        The transformation matrix, T, relating the initial superlattice
        vectors, S, to the prim lattice vectors, L, according to
        ``S = L @ T``, where S and L are shape=(3,3)  matrices with lattice
        vectors as columns.
    dirs: str = "abc"
        The directions along which the initial supercell can be expanded ("a"
        corresponds to first supercell lattice vector, "b" the seoncd, and "c" the
        third).
    min_volume: int = 1000,
        The minimum volume of the results supercell, as integer multiples of the
        prim unit cell volume.

    Returns
    -------
    T_final: array_like, shape=(3,3), dtype=int
        A transformation matrix, ``T_final = T @ M``, scaling T by M, where M is a
        diagonal matrix, such that the volume of `T_final` is greater than or equal to
        `min_volume`.

    """
    T = transformation_matrix_to_super
    M = np.eye(3, dtype=int)
    while int(round(np.linalg.det(T @ M))) < min_volume:
        if "a" in dirs:
            M[0, 0] += 1
        if "b" in dirs:
            M[1, 1] += 1
        if "c" in dirs:
            M[2, 2] += 1
    return T @ M


def min_potential_configs(
    configurations: Union[
        Iterable[casmconfig.Configuration], casmconfig.ConfigurationSet
    ],
    calculator: MonteCalculator,
    conditions: monte.ValueMap,
    transformation_matrix_to_super: Optional[np.ndarray] = None,
    tol: float = casmglobal.TOL,
):
    """Find minimum potential configurations

    Notes
    -----

    - If there is no result (no input configurations, or no input configurations
      that fit the requested supercell), an exception is raised.
    - If there is >1 approximately equivalent minimum potential configuration, the
      first one encountered is returned.

    Parameters
    ----------
    configurations: Union[Iterable[libcasm.configuration.Configuration], \
    libcasm.configuration.ConfigurationSet]
        The candidate configurations. Must be a
        :class:`~libcasm.configuration.ConfigurationSet` or an iterable of
        :class:`~libcasm.configuration.Configuration`.
    calculator: libcasm.clexmonte.MonteCalculator
        A Monte Carlo calculator, to provide access to the potential calculator.
    conditions: Union[dict, libcasm.monte.ValueMap]
        Thermodynamics conditions to calculate the potential at.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        If provided, only configurations that tile the corresponding supercell will
        be considered.
    tol: float = casmglobal.TOL
        Tolerance for identifying configurations with approximately equal potential.

    Returns
    -------
    (configurations, ids, values):

        configurations: list[libcasm.configuration.Configuration]
            The minimum potential configurations.

        ids: list[Union[int, str]]
            If `configurations` is a :class:`~libcasm.configuration.ConfigurationSet`,
            then `ids` is a list of `configuration_name`. If `configurations` is an
            iterable of :class:`~libcasm.configuration.Configuration`, then it is a
            list of indices into the sequence.

        values: list[float]
            The potential per unit cell for the minimum potential configurations.

    """

    # Where we'll store the keep the results as we find them
    min_potential = None
    min_config = []
    min_config_id = []
    min_potential_values = []

    # The MonteCarloState will be set with a configuration when it is checked in
    # order to calculate the potential
    mc_state = None

    # We'll need the prim
    prim = calculator.system.prim

    # If transformation_matrix_to_super: config must tile into mc_big_supercell, but
    # we'll still calculate the potential in config's supercell
    must_tile = False
    if transformation_matrix_to_super is not None:
        mc_big_supercell = casmconfig.Supercell(
            prim=prim,
            transformation_matrix_to_super=transformation_matrix_to_super,
        )
        must_tile = True

    # Loop over candidates
    for i_config, element in enumerate(configurations):
        # 1) check element type, set id and config
        if isinstance(element, casmconfig.Configuration):
            id = i_config
            config = element
        elif isinstance(element, casmconfig.ConfigurationRecord):
            id = element.configuration_name
            config = element.configuration
        else:
            raise Exception(
                "Error in min_potential_config: `configurations` must be an "
                "Iterable[Configuration] or a ConfigurationSet"
            )

        # 2) First time through, set mc_state
        if mc_state is None:
            mc_state = MonteCarloState(
                configuration=config,
                conditions=conditions,
            )
            calculator.set_state_and_potential(state=mc_state)

        # 3) If applicable, skip configs that don't tile into mc_big_supercell
        if must_tile:
            check = mc_big_supercell.superlattice.is_equivalent_superlattice_of(
                config.supercell.superlattice,
                prim.factor_group.elements,
            )
            does_tile, T, fg_index = check
            if not does_tile:
                continue

        # 4) Set the state's configuration
        if config.supercell == mc_state.configuration.supercell:
            # if current config's supercell is the same as the last
            # one we checked, then we don't need to reset the potential calculator,
            # we just copy config DoF values
            mc_state.configuration.dof_values.set(config.dof_values)
        else:
            # otherwise, we need to set the state's configuration
            # and then reset the potential calculator
            mc_state.configuration = config
            calculator.set_state_and_potential(state=mc_state)

        value = calculator.potential.per_unitcell()
        if min_potential is None or value < min_potential - tol:
            min_potential = value
            min_config = [config]
            min_config_id = [id]
            min_potential_values = [value]
        elif value < min_potential + tol:
            min_config.append(config)
            min_config_id.append(id)
            min_potential_values.append(value)

    if len(min_config) == 0:
        raise Exception("Error in min_potential_config: no results")

    return (min_config, min_config_id, min_potential_values)


def make_initial_state(
    calculator: MonteCalculator,
    conditions: Union[dict, monte.ValueMap],
    dirs: str = "abc",
    min_volume: int = 1000,
    transformation_matrix_to_super: Optional[np.ndarray] = None,
    motif: Optional[dict] = None,
    configurations: Union[
        Iterable[casmconfig.Configuration], casmconfig.ConfigurationSet, None
    ] = None,
    tol: float = casmglobal.TOL,
):
    """Determine an appropriate initial state for Monte Carlo calculations

    Parameters
    ----------
    calculator: libcasm.clexmonte.MonteCalculator
        A Monte Carlo calculator, to provide access to the potential calculator.
    conditions: Union[dict, libcasm.monte.ValueMap]
        Thermodynamics conditions to calculate the potential at.
    dirs: str = "abc"
        The directions along which the initial supercell can be expanded ("a"
        corresponds to first supercell lattice vector, "b" the seoncd, and "c" the
        third).
    min_volume: int = 1000,
        The minimum volume of the results supercell, as integer multiples of the
        prim unit cell volume.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        If provided, force the supercell of the result to be a supercell of
        `transformation_matrix_to_supercell`, and only consider motif that tile into
        `transformation_matrix_to_supercell`.
    motif: Optional[libcasm.configuration.Configuration] = None
        If not None, use the provided motif configuration as the initial state rather
        than finding the minimum potential configuration.
    configurations: Union[Iterable[libcasm.configuration.Configuration], \
    libcasm.configuration.ConfigurationSet, None] = None
        The candidate motif configurations. Must be a
        :class:`~libcasm.configuration.ConfigurationSet`, an iterable of
        :class:`~libcasm.configuration.Configuration`, or None. If None, the
        default configuration is used.
    tol: float :data:`~libcasm.casmglobal.TOL`
        Tolerance for identifying configurations with approximately equal potential.

    Returns
    -------
    (initial_state, id):

        initial_state: clexmonte.MonteCarloState
            Initial Monte Carlo state according to the specified parameters.
        id: Union[str, int]
            ID of the configuration chosen to fill `initial_state`. May be ``"motif"``,
            if a motif was provided, ``"default"`` if no motif or configurations were
            provided, a `configuration_name` string if a ConfigurationSet was provided,
            or an integer index into the sequence if an iterable of Configuration was
            provided.
    """
    if motif is None:
        add_default = False
        if configurations is None:
            add_default = True
            supercell = casmconfig.Supercell(
                prim=calculator.system.prim,
                transformation_matrix_to_super=np.eye(3, dtype="int"),
            )
            default_configuration = casmconfig.Configuration(
                supercell=supercell,
            )
            configurations = [default_configuration]

        configs, ids, values = min_potential_configs(
            configurations=configurations,
            calculator=calculator,
            conditions=conditions,
            transformation_matrix_to_super=transformation_matrix_to_super,
            tol=tol,
        )

        # remove equivalent configurations
        if len(configs) > 1:
            _prim_canonical_configs = []
            _configs = []
            _ids = []
            _values = []
            for i in range(len(configs)):
                x = casmconfig.make_primitive_configuration(configs[i])
                x = casmconfig.make_canonical_configuration(
                    x, in_canonical_supercell=True
                )
                if x not in _prim_canonical_configs:
                    _prim_canonical_configs.append(x)
                    _configs.append(configs[i])
                    _ids.append(ids[i])
                    _values.append(values[i])
            configs = _configs
            ids = _ids
            values = _values

        # warn if >1 result found
        if len(configs) > 1:
            print(
                "Warning: auto_configuration found >1 configuration with potential "
                "approximately equal to the minimum"
            )
            for i in range(len(configs)):
                print(
                    f"- id: {ids[i]}, "
                    f"potential: {values[i]}, "
                    f"occ: {configs[i].occupation}"
                )
            print(f"Using: id={ids[0]}")
        motif = configs[0]

        if add_default:
            id = "default"
        else:
            id = ids[0]
    else:
        id = "motif"

    if transformation_matrix_to_super is None:
        T_init = motif.supercell.transformation_matrix_to_super
    else:
        T_init = transformation_matrix_to_super

    T = scale_supercell(
        transformation_matrix_to_super=T_init,
        dirs=dirs,
        min_volume=min_volume,
    )
    supercell = casmconfig.Supercell(
        prim=motif.supercell.prim,
        transformation_matrix_to_super=T,
    )
    check = supercell.superlattice.is_equivalent_superlattice_of(
        motif.supercell.superlattice,
        motif.supercell.prim.factor_group.elements,
    )
    does_tile, T, fg_index = check
    if not does_tile:
        raise Exception(
            "Error in make_initial_state: Failed to tile motif into supercell"
        )
    config = casmconfig.copy_transformed_configuration(
        prim_factor_group_index=fg_index,
        translation=[0, 0, 0],
        motif=motif,
        supercell=supercell,
    )
    return (
        MonteCarloState(
            configuration=config,
            conditions=conditions,
        ),
        motif,
        id,
    )
