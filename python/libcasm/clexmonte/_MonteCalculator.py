import textwrap

from ._clexmonte_monte_calculator import (
    MonteCalculator,
)


def print_selected_event_functions(calculator: MonteCalculator):
    """Print a summary of the selected event functions in a MonteCalculator

    Parameters
    ----------
    calculator : MonteCalculator
        The MonteCalculator to summarize the selected event functions of

    """

    all_functions = calculator.selected_event_functions

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
                fill("  Requires event state = " + str(function.requires_event_state))
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
                    print(fill("  Component names = " + str(function.component_names)))
            if hasattr(function, "partition_names"):
                print(fill("  Partition names = " + str(function.partition_names)))
            if hasattr(function, "value_labels"):
                value_labels = function.value_labels()
                if value_labels is not None:
                    labels = [x[1] for x in value_labels]
                    print(fill("  Value labels = " + str(labels)))
            print(
                fill("  Requires event state = " + str(function.requires_event_state))
            )
            if hasattr(function, "is_log"):
                if function.is_log:
                    print(fill("  Is log = " + str(function.is_log)))
            if hasattr(function, "initial_begin"):
                print(fill("  Default initial bin = " + str(function.initial_begin)))
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
