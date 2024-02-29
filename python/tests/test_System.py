import libcasm.xtal as xtal
import libcasm.configuration as casmconfig
from libcasm.clexmonte import (
    System,
)
from libcasm.clexulator import (
    PrimNeighborList,
)
from libcasm.composition import (
    CompositionCalculator,
    CompositionConverter,
)


def test_System_constructor_1(
        FCCBinaryVacancy_xtal_prim,
        FCCBinaryVacancy_CompositionConverter,
):
    system = System(
        xtal_prim=FCCBinaryVacancy_xtal_prim,
        composition_converter=FCCBinaryVacancy_CompositionConverter,
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_System_from_dict_1():
    data = {
        "prim": {
            "basis": [
                {
                    "coordinate": [0.0, 0.0, 0.0],
                    "occupants": ["A", "B", "Va"]
                }
            ],
            "coordinate_mode": "Fractional",
            "lattice_vectors": [
                [0.0, 2.0, 2.0],
                [2.0, 0.0, 2.0],
                [2.0, 2.0, 0.0]
            ],
            "title": "FCC_binary_vacancy"
        },
        "composition_axes": {
            "a": [0.000000000000, 1.000000000000, 0.000000000000],
            "b": [0.000000000000, 0.000000000000, 1.000000000000],
            "components": ["A", "B", "Va"],
            "independent_compositions": 2,
            "mol_formula": "A(1-a-b)B(a)Va(b)",
            "origin": [1.000000000000, 0.000000000000, 0.000000000000],
            "param_formula": "a(0.333333-0.333333A+0.666667B-0.333333Va)b(0.333333-0.333333A-0.333333B+0.666667Va)"
        },
    }
    system = System.from_dict(data)
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)
