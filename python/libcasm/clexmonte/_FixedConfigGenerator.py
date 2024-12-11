import copy
from typing import Optional

import numpy as np

import libcasm.monte as monte
from libcasm.configuration import (
    Configuration,
    Supercell,
    SupercellSet,
    copy_transformed_configuration,
)

from ._RunData import (
    RunData,
)
from .parsing import (
    optional_from_dict,
    to_dict,
)


class FixedConfigGenerator:
    """A `ConfigGenerator` for state generation - always returns the same configuration

    Notes
    -----

    - Returns the same configuration no matter what the current conditions and
      completed runs are.

    """

    def __init__(
        self,
        configuration: Optional[Configuration] = None,
        supercell: Optional[Supercell] = None,
        motif: Optional[Configuration] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        configuration:  Optional[libcasm.configuration.Configuration]=None
            The configuration to generate, directly. If provided, other
            parameters are ignored.
        supercell: Optional[libcasm.configuration.Supercell]=None
            The Monte Carlo supercell if `configuration` is not provided.
        motif: Optional[libcasm.configuration.Configuration]=None
            Initial configuration, which will be copied and tiled into the
            Monte Carlo supercell. If a perfect tiling can be made by
            applying factor group operations, a note is printed indicating
            which operation is applied. A warning is printed if there is
            no perfect tiling and the `motif` is used without reorientation
            to fill the supercell imperfectly. If `supercell` is given but
            no `motif` is provided, the default configuration is used.

        """
        self._configuration = None
        self._supercell = None
        self._motif = None
        self._default_motif = None

        if configuration is not None:
            self._configuration = copy.copy(configuration)
        elif supercell is not None:
            self._supercell = supercell
            if motif is None:
                self._motif = Configuration(self._supercell)
                self._default_motif = True
                fg_index = 0
            else:
                self._motif = motif
                self._default_motif = False
                factor_group = supercell.prim.factor_group
                (
                    is_equivalent,
                    T,
                    fg_index,
                ) = supercell.superlattice.is_equivalent_superlattice_of(
                    self._motif.superlattice, factor_group.elements
                )
                if not is_equivalent:
                    print(
                        "Warning: `motif` cannot tile `supercell`. "
                        "Will fill imperfectly."
                    )
                elif fg_index != 0:
                    print(
                        f"Note: `motif` fills `supercell` after applying factor "
                        f"group operation {fg_index} (indexing from 0)."
                    )
            self._configuration = copy_transformed_configuration(
                prim_factor_group_index=fg_index,
                translation=[0, 0, 0],
                motif=self._motif,
                supercell=self._supercell,
            )
        else:
            raise Exception(
                "Error constructing FixedConfigGenerator: "
                "One of `configuration` or `supercell` must be provided."
            )

    def __call__(
        self,
        conditions: monte.ValueMap,
        completed_runs: list[RunData],
    ) -> Configuration:
        """Always returns the same configuration

        Parameters
        ----------
        conditions: libcasm.monte.ValueMap
            The thermodynamic conditions at which the next state will be run
        completed_runs: list[RunData]
            Ignored, but included for compatibility.

        Returns
        -------
        configuration:  libcasm.configuration.Configuration
            The configuration provided to the constructor.
        """
        return copy.copy(self._configuration)

    @staticmethod
    def methodname() -> str:
        """Returns the method name "fixed" """
        return "fixed"

    def to_dict(
        self,
    ):
        """Convert FixedConfigGenerator to a Python dict

        Returns
        -------
        data: dict
            The FixedConfigGenerator as a Python dict
        """
        if self._motif is not None:
            data = {}
            to_dict(
                self._supercell.transformation_matrix_to_super.tolist(),
                data,
                "transformation_matrix_to_supercell",
            )
            if self._default_motif is False:
                to_dict(self._motif, data, "motif")
            return data
        else:
            data = {}
            to_dict(self._configuration, data, "configuration")
            return data

    @staticmethod
    def from_dict(
        data: dict,
        supercells: SupercellSet,
    ):
        """Construct FixedConfigGenerator from a Python dict

        Returns
        -------
        config_generator: FixedConfigGenerator
            The FixedConfigGenerator
        """
        configuration = optional_from_dict(
            Configuration, data, "configuration", supercells=supercells
        )
        supercell = None
        transformation_matrix_to_super = data.get("transformation_matrix_to_supercell")
        if transformation_matrix_to_super is not None:
            T = np.array(transformation_matrix_to_super, dtype="int")
            supercell = supercells.add_by_transformation_matrix_to_super(T).supercell
        motif = optional_from_dict(Configuration, data, "motif", supercells=supercells)

        return FixedConfigGenerator(
            configuration=configuration, supercell=supercell, motif=motif
        )
