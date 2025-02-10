"""Find minimum potential configurations"""

from collections.abc import Iterable
from typing import Optional, Union

import numpy as np

import libcasm.casmglobal as casmglobal
import libcasm.configuration as casmconfig
import libcasm.monte as monte
import libcasm.xtal as xtal

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


def transform_configuration(
    prim_factor_group_index: int,
    config: casmconfig.Configuration,
):
    prim = config.supercell.prim
    symop = prim.factor_group.elements[prim_factor_group_index]
    P = prim.xtal_prim.lattice()
    S1 = config.supercell.superlattice
    S2 = xtal.make_canonical_lattice(symop * S1)
    is_superlattice_of, T2 = S2.is_superlattice_of(P)
    if not is_superlattice_of:
        raise Exception(
            "Error in transform_configuration: "
            "construction of transformed supercell failed"
        )
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.rint(T2).astype(int),
    )
    return casmconfig.copy_transformed_configuration(
        prim_factor_group_index=prim_factor_group_index,
        translation=[0, 0, 0],
        motif=config,
        supercell=supercell,
        origin=[0, 0, 0],
    )


class FindMinPotentialConfigs:
    """Find minimum potential configurations"""

    def __init__(
        self,
        calculator: MonteCalculator,
        conditions: monte.ValueMap,
        transformation_matrix_to_super: Optional[np.ndarray] = None,
        isotropic_potential: bool = True,
        tol: float = casmglobal.TOL,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        calculator: libcasm.clexmonte.MonteCalculator
            A Monte Carlo calculator, to provide access to the potential calculator.
        conditions: Union[dict, libcasm.monte.ValueMap]
            Thermodynamics conditions to calculate the potential at.
        transformation_matrix_to_super: Optional[np.ndarray] = None
            If provided, only configurations that tile the corresponding supercell will
            be considered.
        isotropic_potential: bool = True
            If True, search for minimum potential configurations assuming that
            orientation does not matter. If False, apply the prim factor group to
            check all possible orientations.
        tol: float = casmglobal.TOL
            Tolerance for identifying configurations with approximately equal potential.

        """

        self.calculator = calculator
        self.prim = self.calculator.system.prim
        self.conditions = conditions
        self.transformation_matrix_to_super = transformation_matrix_to_super
        self.isotropic_potential = isotropic_potential
        self.tol = tol

        # If transformation_matrix_to_super: config must tile into supercell_1, but
        # we'll still calculate the potential in config's supercell
        self.must_tile = False
        self.supercell_1 = None
        if self.transformation_matrix_to_super is not None:
            self.supercell_1 = casmconfig.Supercell(
                prim=self.prim,
                transformation_matrix_to_super=self.transformation_matrix_to_super,
            )
            self.must_tile = True

        self.reset()

    def reset(self):
        """Reset stored results"""
        # Where we'll store the results as we find them
        self.min_potential = None
        self.min_config = []
        self.min_config_id = []
        self.min_potential_values = []

        # The MonteCarloState will be set with a configuration when it is checked in
        # order to calculate the potential
        self._mc_state = None

    def check(
        self,
        configurations: Union[
            Iterable[casmconfig.Configuration], casmconfig.ConfigurationSet
        ],
    ):
        """Check for minimum potential configurations

        Notes
        -----

        - If there is no result (empty input configurations, or no input configurations
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

        Returns
        -------
        (configurations, ids, values):

            configurations: list[libcasm.configuration.Configuration]
                The minimum potential configurations.

            ids: list[Union[int, str]]
                If `configurations` is a
                :class:`~libcasm.configuration.ConfigurationSet`, then `ids` is a list
                of `configuration_name`. If `configurations` is an iterable of
                :class:`~libcasm.configuration.Configuration`, then it is a list of
                indices into the sequence.

            values: list[float]
                The potential per unit cell for the minimum potential configurations.

        """
        # Loop over candidates
        for i_config, element in enumerate(configurations):
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
            self._check_config(config, id)

        if len(self.min_config) == 0:
            raise Exception("Error in FindMinPotentialConfigs.check: no results")

        return (self.min_config, self.min_config_id, self.min_potential_values)

    def _check_config(
        self,
        config: casmconfig.Configuration,
        id: str,
    ):
        if self.isotropic_potential:
            if self.must_tile:
                check = self.supercell_1.superlattice.is_equivalent_superlattice_of(
                    config.supercell.superlattice,
                    self.prim.factor_group.elements,
                )
                does_tile, T, prim_factor_group_index = check
                if not does_tile:
                    return
                config_orientation = transform_configuration(
                    prim_factor_group_index, config
                )
                self._check_single_orientation(config_orientation, id)
            else:
                self._check_single_orientation(config, id)
        else:
            for prim_factor_group_index in range(len(self.prim.factor_group.elements)):
                config_orientation = transform_configuration(
                    prim_factor_group_index, config
                )
                self._check_single_orientation(config_orientation, id)

    def _check_single_orientation(
        self,
        config: casmconfig.Configuration,
        id: str,
    ):
        # If not set, set _mc_state and potential calculator
        if self._mc_state is None:
            self._mc_state = MonteCarloState(
                configuration=config,
                conditions=self.conditions,
            )
            self.calculator.set_state_and_potential(state=self._mc_state)

        # If applicable, skip configs that don't tile into mc_big_supercell
        if self.must_tile:
            is_superlattice_of, T = self.supercell_1.superlattice.is_superlattice_of(
                config.supercell.superlattice
            )
            if not is_superlattice_of:
                return

        # Set the state's configuration
        if config.supercell == self._mc_state.configuration.supercell:
            # if current config's supercell is the same as the last
            # one we checked, then we don't need to reset the potential calculator,
            # we just copy config DoF values
            self._mc_state.configuration.dof_values.set(config.dof_values)
        else:
            # otherwise, we need to set the state's configuration
            # and then reset the potential calculator
            self._mc_state.configuration = config
            self.calculator.set_state_and_potential(state=self._mc_state)

        # Calculate the potential value
        value = self.calculator.potential.per_unitcell()
        if self.min_potential is None or value < self.min_potential - self.tol:
            self.min_potential = value
            self.min_config = [config]
            self.min_config_id = [id]
            self.min_potential_values = [value]
        elif value < self.min_potential + self.tol:
            self.min_config.append(config)
            self.min_config_id.append(id)
            self.min_potential_values.append(value)


def make_initial_state(
    calculator: MonteCalculator,
    conditions: Union[dict, monte.ValueMap],
    dirs: str = "abc",
    min_volume: int = 1000,
    transformation_matrix_to_super: Optional[np.ndarray] = None,
    motifs: Union[
        Iterable[casmconfig.Configuration], casmconfig.ConfigurationSet, None
    ] = None,
    isotropic_potential: bool = True,
    tol: float = casmglobal.TOL,
):
    """Determine an appropriate initial state for Monte Carlo calculations

    Notes
    -----

    The method here works according to the following steps:

    1. Find a motif configuration:

       - Use the calculator's potential to find the motif configuration at the
         given `conditions`, conditioned on the following
         restrictions:

         - If `motifs` is given, the search is made over orientations of the provided
           configurations; else only the default configuration in the primitive unit
           cell is considered.
         - If `transformation_matrix_to_super` is given, the chosen motif must tile
           into the supercell specified by `transformation_matrix_to_super`.
         - The prim factor group is used to generate orientations of the
           configurations under consideration. If `isotropic_potential` is True, only
           the first orientation which meets the tiling criteria is considered. If
           `isotropic_potential` is False, all orientations which meet the tiling
           criteria are checked individually.
         - If multiple candidate motifs have approximately the same minimum potential,
           then:

           - Equivalent configurations are removed, as determined by having the same
             primitive, canonical form;
           - If there are still multiple candidate motifs, then a warning is printed
             listing all the equivalents and the first one found during the search
             is chosen as the motif.

    2. Find the initial state's supercell

       - First, a minimum supercell is constructed:

         - If `transformation_matrix_to_super` is provided, it is used to construct the
           minimum supercell;
         - Else, the supercell of the motif found at the previous step is used as the
           minimum supercell.

       - Then the minimum supercell is scaled using :func:`scale_supercell`, subject to
         `dirs` and `min_volume`, to construct the initial state's supercell

    3. Fill the supercell with the motif configuration

       - The chosen motif is tiled into the initial state's supercell to create the
         initial state.


    Parameters
    ----------
    calculator: libcasm.clexmonte.MonteCalculator
        A Monte Carlo calculator, to provide access to the potential calculator.
    conditions: Union[dict, libcasm.monte.ValueMap]
        Thermodynamics conditions to calculate the potential at.
    dirs: str = "abc"
        The directions along which the initial supercell can be expanded ("a"
        corresponds to first supercell lattice vector, "b" the second, and "c" the
        third).
    min_volume: int = 1000,
        The minimum volume of the results supercell, as integer multiples of the
        prim unit cell volume.
    transformation_matrix_to_super: Optional[np.ndarray] = None
        If provided, force the supercell of the result to be a supercell of
        `transformation_matrix_to_supercell`, and only consider motif that tile into
        `transformation_matrix_to_supercell`.
    motifs: Union[Iterable[libcasm.configuration.Configuration], \
    libcasm.configuration.ConfigurationSet, None] = None
        The candidate motif configurations. Must be a
        :class:`~libcasm.configuration.ConfigurationSet`, an iterable of
        :class:`~libcasm.configuration.Configuration`, or None. If None, the
        default configuration is used.
    isotropic_potential: bool = True
        If True, search for minimum potential configurations assuming that
        orientation does not matter. If False, apply the prim factor group to
        check all possible orientations.
    tol: float :data:`~libcasm.casmglobal.TOL`
        Tolerance for identifying configurations with approximately equal potential.

    Returns
    -------
    initial_state: libcasm.clexmonte.MonteCarloState
        Initial Monte Carlo state according to the specified parameters.
    motif: libcasm.configuration.Configuration
        The motif chosen to fill `initial_state`.
    id: Union[str, int]
        ID of the motif configuration chosen to fill `initial_state`. May be
        ``"default"`` if no motif configurations were provided, a
        `configuration_name` string if a ConfigurationSet was provided,
        or an integer index into the sequence if an iterable of Configuration was
        provided.
    """

    add_default = False
    if motifs is None:
        add_default = True
        supercell = casmconfig.Supercell(
            prim=calculator.system.prim,
            transformation_matrix_to_super=np.eye(3, dtype="int"),
        )
        default_configuration = casmconfig.Configuration(
            supercell=supercell,
        )
        motifs = [default_configuration]

    finder = FindMinPotentialConfigs(
        calculator=calculator,
        conditions=conditions,
        transformation_matrix_to_super=transformation_matrix_to_super,
        isotropic_potential=isotropic_potential,
        tol=tol,
    )

    configs, ids, values = finder.check(configurations=motifs)

    # remove equivalent configurations
    if len(configs) > 1:
        _prim_canonical_configs = []
        _configs = []
        _ids = []
        _values = []
        for i in range(len(configs)):
            x = casmconfig.make_primitive_configuration(configs[i])
            x = casmconfig.make_canonical_configuration(x, in_canonical_supercell=True)
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
            "Warning: auto_configuration found >1 symmetrically distinct "
            "configuration with potential approximately equal to the minimum."
        )
        report = []
        for i in range(len(configs)):
            report.append(
                {
                    "id": ids[i],
                    "potential": values[i],
                    "config": configs[i].to_dict(),
                }
            )
        print("Minimum potential motifs:")
        print(xtal.pretty_json(report))
        print("Using:")
        print(xtal.pretty_json(report[0]))
    motif = configs[0]

    if add_default:
        id = "default"
    else:
        id = ids[0]

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
    is_superlattice_of, T = supercell.superlattice.is_superlattice_of(
        motif.supercell.superlattice
    )
    if not is_superlattice_of:
        raise Exception(
            "Error in make_initial_state: Failed to tile motif into supercell"
        )
    config = casmconfig.copy_configuration(
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


def make_canonical_initial_state(
    calculator: MonteCalculator,
    conditions: Union[dict, monte.ValueMap],
    dirs: str = "abc",
    min_volume: int = 1000,
    transformation_matrix_to_super: Optional[np.ndarray] = None,
    motif: Optional[dict] = None,
    isotropic_potential: bool = True,
    tol: float = casmglobal.TOL,
):
    """Determine an appropriate initial state for canonical Monte Carlo calculations

    Notes
    -----

    The method here works according to the following steps:

    1. With the `motif` provided, or using the default configuration in the primitive
       unit cell, use `make_initial_state` to construct an initial Monte Carlo state,
       optimally oriented in the supercell, satisfying `transformation_matrix_to_super`,
       `dirs`, and `min_volume`.

    2. If either `mol_composition` or `param_composition` are included in `conditions`,
       modify the configuration of the state generated by the previous step to
       enforce the specified composition as closely as possible. Modification is
       performed using the calculator's `"enforce.composition"` modifying function. For
       standard MonteCalculator implementations, the standard `"enforce.composition"`
       function determines which type(s) of occupant swap could be performed to adjust
       the composition most rapidly towards the target value and uses
       `calculator.engine` to randomly choose which swap to make.

    3. Return the Monte Carlo state with configuration generated by the previous step
       and conditions generated by taking the input `conditions` and modifying the
       composition to exactly match the composition of the configuration, using the
       calculator's `"match.composition"` modifying function.

    Parameters
    ----------
    calculator: libcasm.clexmonte.MonteCalculator
        A Monte Carlo calculator, to provide access to the potential calculator. The
        calculator must have state modifying functions named `"enforce.composition"`
        (which modifies the configuration to match the composition conditions) and
        `"match.composition"` (which modifies the composition conditions to match the
        configuration).
    conditions: Union[dict, libcasm.monte.ValueMap]
        Thermodynamics conditions for the initial state. If neither `mol_composition`
        nor `param_composition` are included, then they are calculated from the
        input motif.
    dirs: str = "abc"
        The directions along which the initial supercell can be expanded ("a"
        corresponds to first supercell lattice vector, "b" the second, and "c" the
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
    isotropic_potential: bool = True
        If True, search for minimum potential configurations assuming that
        orientation does not matter. If False, apply the prim factor group to
        check all possible orientations.
    tol: float :data:`~libcasm.casmglobal.TOL`
        Tolerance for identifying configurations with approximately equal potential.

    Returns
    -------
    initial_state: clexmonte.MonteCarloState
        Initial Monte Carlo state according to the specified parameters.
    motif: libcasm.configuration.Configuration
        The motif chosen to fill `initial_state`.
    """
    match_composition_f = calculator.modifying_functions["match.composition"]
    enforce_composition_f = calculator.modifying_functions["enforce.composition"]

    if motif is None:
        supercell = casmconfig.Supercell(
            prim=calculator.system.prim,
            transformation_matrix_to_super=np.eye(3, dtype="int"),
        )
        motif = casmconfig.Configuration(
            supercell=supercell,
        )

    # mc_state_1:
    # - input motif,
    # - conditions w/ composition matching input motif
    mc_state_1 = MonteCarloState(
        configuration=motif,
        conditions=conditions,
    )
    match_composition_f(mc_state_1)
    if "mol_composition" not in conditions and "param_composition" not in conditions:
        conditions = mc_state_1.conditions

    # mc_state_2:
    # - best oriented motif,
    # - conditions w/ composition matching input motif
    mc_state_2, chosen_motif, id = make_initial_state(
        calculator=calculator,
        conditions=mc_state_1.conditions,
        dirs=dirs,
        min_volume=min_volume,
        transformation_matrix_to_super=transformation_matrix_to_super,
        motifs=[motif],
        isotropic_potential=isotropic_potential,
        tol=tol,
    )

    # mc_state_3:
    # - best oriented motif modified to enforce composition of input conditions,
    # - conditions w/ composition set to exactly match the enforced composition
    mc_state_3 = MonteCarloState(
        configuration=mc_state_2.configuration,
        conditions=conditions,
    )
    enforce_composition_f(mc_state_3)
    match_composition_f(mc_state_3)

    return (mc_state_3, chosen_motif)
