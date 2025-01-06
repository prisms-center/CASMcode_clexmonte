import textwrap
from typing import Optional

from libcasm.local_configuration import (
    LocalConfiguration,
)

from ._clexmonte_monte_calculator import (
    MonteCalculatorCore,
)
from ._clexmonte_system import (
    System,
)
from ._system_methods import (
    make_system_event_info,
)


class MonteCalculator(MonteCalculatorCore):
    def __init__(
        self,
        method: str,
        system: System,
        params: Optional[dict] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        method : str
            Monte Carlo method name. The options are:

            - "semigrand_canonical": Metropolis algorithm in the semi-grand
              canonical ensemble. Input states require `"temperature"` and
              `"param_chem_pot"` conditions.
            - "canonical": Metropolis algorithm in the canonical ensemble.
              Input states require `"temperature"` and one of
              `"param_composition"` or `"mol_composition"` conditions.
            - "kinetic": Kinetic Monte Carlo method. Input states require
              `"temperature"` and one of `"param_composition"` or
              `"mol_composition"` conditions.

        system : libcasm.clexmonte.System
            Cluster expansion model system data. The required data depends on
            the calculation method. See links under `method` for what system
            data is required for each method.

        params: Optional[dict] = None
            Monte Carlo calculation method parameters. Expected values
            depends on the calculation method.

        """
        super().__init__(
            method=method,
            system=system,
            params=params,
        )

        self.event_info = make_system_event_info(self.system)
        """dict[str, libcasm.local_configuration.OccEventSymInfo]: A dictionary of 
        event info used to construct and transform 
        :class:`~libcasm.local_configuration.LocalConfiguration`."""

    def make_local_configuration(self):
        """Make a :class:`~libcasm.local_configuration.LocalConfiguration` from the
        current state and current selected event.

        Returns
        -------
        local_configuration : libcasm.local_configuration.LocalConfiguration
            The local configuration.
        """
        selected_event = self.selected_event
        event_type_name = selected_event.prim_event_data.event_type_name
        unitcell_index = selected_event.event_data.unitcell_index
        equivalent_index = selected_event.prim_event_data.equivalent_index
        return LocalConfiguration(
            configuration=self.state_data.state.configuration,
            pos=(unitcell_index, equivalent_index),
            event_info=self.event_info[event_type_name],
        )

    def print_selected_event_functions(self):
        """Print a summary of the selected event functions in a MonteCalculator"""

        all_functions = self.selected_event_functions

        def fill(text):
            return textwrap.fill(
                text,
                width=80,
                initial_indent="",
                subsequent_indent="    ",
            )

        def print_generic_functions(functions):
            for key, function in functions.items():
                print(key + ":")
                print(fill("  Description = " + function.description))
                print(
                    fill(
                        "  Requires event state = " + str(function.requires_event_state)
                    )
                )
                print(fill("  Default order = " + str(function.order)))
                print()

        def print_functions(functions):
            for key, function in functions.items():
                print(key + ":")
                print(fill("  Description = " + function.description))
                if hasattr(function, "shape"):
                    if len(function.shape) == 0:
                        print(fill("  Shape = [] (Scalar)"))
                    else:
                        print(fill("  Shape = " + str(function.shape)))
                        print(
                            fill("  Component names = " + str(function.component_names))
                        )
                if hasattr(function, "partition_names"):
                    print(fill("  Partition names = " + str(function.partition_names)))
                if hasattr(function, "value_labels"):
                    value_labels = function.value_labels()
                    if value_labels is not None:
                        labels = [x[1] for x in value_labels]
                        print(fill("  Value labels = " + str(labels)))
                print(
                    fill(
                        "  Requires event state = " + str(function.requires_event_state)
                    )
                )
                if hasattr(function, "is_log"):
                    if function.is_log:
                        print(fill("  Is log = " + str(function.is_log)))
                if hasattr(function, "initial_begin"):
                    print(
                        fill("  Default initial bin = " + str(function.initial_begin))
                    )
                    print(fill("  Default bin width = " + str(function.bin_width)))
                print(fill("  Default max size = " + str(function.max_size)))
                if hasattr(function, "tol"):
                    print(fill("  Default tol = " + str(function.tol)))
                print()

        functions = all_functions.generic_functions
        if len(functions):
            print("Selected event functions:\n")
            print_generic_functions(functions)

        int_functions = all_functions.discrete_vector_int_functions
        float_functions = all_functions.discrete_vector_float_functions
        continuous_1d_functions = all_functions.continuous_1d_functions

        if len(int_functions) + len(float_functions) + len(continuous_1d_functions):
            print("Selected event data functions:\n")
            print_functions(int_functions)
            print_functions(float_functions)
            print_functions(continuous_1d_functions)
