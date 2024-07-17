from typing import Optional

import numpy as np

from libcasm.configuration import (
    SupercellSet,
)
from libcasm.monte import ValueMap

from ._clexmonte_state import MonteCarloState
from .parsing import (
    optional_from_dict,
    required_from_dict,
    required_int_array_from_dict,
    to_dict,
)


class RunData:
    """Data structure to store data for each run in a series

    RunData is used by :class:`~libcasm.clexmonte.IncrementalConditionsStateGenerator`
    for:

    - storing a summary of the runs in a series,
    - determining the next run in the series, including its initial state in the
      case of dependent runs, and
    - enabling restarts.

    """

    def __init__(
        self,
        initial_state: MonteCarloState,
        final_state: Optional[MonteCarloState] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        initial_state: MonteCarloState
            Initial Monte Carlo state
        final_state: Optional[MonteCarloState] = None
            Final Monte Carlo state
        """
        self.transformation_matrix_to_super = (
            initial_state.configuration.supercell.transformation_matrix_to_super
        )
        """np.ndarray[np.int64]: Monte Carlo supercell
        
        The transformation matrix, T, relating the Monte Carlo supercell lattice
        vectors, S, to the primitive cell lattice vectors, L, according to
        ``S = L @ T``, where S and L are shape=(3,3)  matrices with
        lattice vectors as columns.
        """

        self.n_unitcells = int(
            round(np.linalg.det(self.transformation_matrix_to_super))
        )
        """int: Monte Carlo supercell volume, as integer multiple of the primitive cell
        """

        self.conditions = initial_state.conditions
        """ValueMap: Thermodynamic conditions"""

        self.initial_state = initial_state
        """Optional[MonteCarloState]: Initial Monte Carlo state"""

        self.final_state = final_state
        """Optional[MonteCarloState]: Final Monte Carlo state"""

    def to_dict(
        self,
        write_initial_states: bool = False,
        write_final_states: bool = False,
        write_prim_basis: bool = False,
    ):
        """Convert RunData to a Python dict

        Parameters
        ----------
        write_initial_states: bool = False
            If True, write initial_state.
        write_final_states: bool = False
            If True, write final_state.
        write_prim_basis: bool = False
            If True, write state configuration degrees of freedom (DoF) using prim
            basis. Default is to use the standard basis.

        Returns
        -------
        data: dict
            The RunData as a Python dict
        """
        data = {}
        to_dict(self.conditions, data, "conditions")
        to_dict(
            self.transformation_matrix_to_super.tolist(),
            data,
            "transformation_matrix_to_supercell",
        )
        to_dict(self.n_unitcells, data, "n_unitcells")
        to_dict(self.conditions, data, "conditions")
        if write_initial_states:
            to_dict(
                self.initial_state,
                data,
                "initial_state",
                write_prim_basis=write_prim_basis,
            )
        if write_final_states:
            to_dict(
                self.final_state,
                data,
                "final_state",
                write_prim_basis=write_prim_basis,
            )
        return data

    @staticmethod
    def from_dict(
        data: dict,
        supercells: SupercellSet,
    ):
        """Construct RunData from a Python dict"""
        return RunData(
            transformation_matrix_to_super=required_int_array_from_dict(
                data, "transformation_matrix_to_supercell"
            ),
            n_unitcells=required_from_dict(int, data, "n_unitcells"),
            conditions=required_from_dict(ValueMap, data, "conditions"),
            initial_state=optional_from_dict(
                MonteCarloState, data, "initial_state", supercells=supercells
            ),
            final_state=optional_from_dict(
                MonteCarloState, data, "final_state", supercells=supercells
            ),
        )


class RunDataOutputParams:
    """Parameters controlling saving and write run data for a series of runs

    RunDataOutputParams is used by
    :class:`~libcasm.clexmonte.IncrementalConditionsStateGenerator` to control:

    - storing a summary of the runs in a series,
    - determining the next run in the series, including its initial state in the
      case of dependent runs, and
    - enabling restarts.

    """

    def __init__(
        self,
        do_save_all_initial_states: bool = False,
        do_save_all_final_states: bool = False,
        do_save_last_final_state: bool = True,
        write_initial_states: bool = False,
        write_final_states: bool = False,
        write_prim_basis: bool = False,
        output_dir: Optional[str] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        do_save_all_initial_states: bool = False
            Save all initial states in the state generator's completed runs list
        do_save_all_final_states: bool = False
            Save all final states in the state generator's completed runs list
        do_save_last_final_state: bool = True
            Save the last final state in the state generator's completed runs list (to
            enable restarts), and if `write_final_states` is False, delete previous
            saved final states to reduce memory usage and file size.
        write_initial_states: bool = False
            Write saved initial states to completed_runs.json
        write_final_states: bool = False
            Write saved final states to completed_runs.json
        write_prim_basis: bool = False
            Write saved state degrees of freedom (DoF) values using the prim basis.
            Default (False) is to use the standard basis.
        output_dir: Optional[str] = None
            Location to save completed_runs.json if not None
        """
        self.do_save_all_initial_states = do_save_all_initial_states
        """bool: Save all initial states in the state generator's completed runs list"""

        self.do_save_all_final_states = do_save_all_final_states
        """bool: Save all final states in the state generator's completed runs list"""

        self.do_save_last_final_state = do_save_last_final_state
        """bool: Save the last final state in the state generator's completed runs list 
        (to enable restarts), and if `write_final_states` is False, delete previous
        saved final states to reduce memory usage and file size."""

        self.write_initial_states = write_initial_states
        """bool: Write saved initial states to completed_runs.json"""

        self.write_final_states = write_final_states
        """bool: Write saved final states to completed_runs.json"""

        self.write_prim_basis = write_prim_basis
        """bool: Write saved states degrees of freedom (DoF) values using the prim \
        basis."""

        self.output_dir = output_dir
        """Optional[str]: Location to save completed_runs.json if not None"""

    def to_dict(
        self,
        write_initial_states: bool = False,
        write_final_states: bool = False,
    ):
        """Convert RunDataOutputParams to a Python dict"""
        data = {}
        to_dict(self.do_save_all_initial_states, data, "save_all_initial_states")
        to_dict(self.do_save_all_final_states, data, "save_all_final_states")
        to_dict(self.do_save_last_final_state, data, "save_last_final_state")
        to_dict(self.write_initial_states, data, "write_initial_states")
        to_dict(self.write_final_states, data, "write_final_states")
        to_dict(self.write_prim_basis, data, "write_prim_basis")
        to_dict(self.output_dir, data, "output_dir")
        return data

    @staticmethod
    def from_dict(
        data: dict,
    ):
        """Construct RunDataOutputParams from a Python dict"""
        return RunDataOutputParams(
            do_save_all_initial_states=optional_from_dict(
                bool, data, "save_all_initial_states", default_value=False
            ),
            do_save_all_final_states=optional_from_dict(
                bool, data, "save_all_final_states", default_value=False
            ),
            do_save_last_final_state=optional_from_dict(
                bool, data, "save_last_final_state", default_value=True
            ),
            write_initial_states=optional_from_dict(
                bool, data, "write_initial_states", default_value=False
            ),
            write_final_states=optional_from_dict(
                bool, data, "write_final_states", default_value=False
            ),
            write_prim_basis=optional_from_dict(
                bool, data, "write_prim_basis", default_value=False
            ),
            output_dir=optional_from_dict(str, data, "output_dir", default_value=None),
        )
