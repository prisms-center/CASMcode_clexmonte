import pytest
import libcasm.xtal as xtal
from libcasm.composition import (
    CompositionConverter,
)


@pytest.fixture
def FCCBinaryVacancy_xtal_prim():
    data = {
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
    }
    return xtal.Prim.from_dict(data)


@pytest.fixture
def FCCBinaryVacancy_CompositionConverter():
    data = {
        "a": [0.000000000000, 1.000000000000, 0.000000000000],
        "b": [0.000000000000, 0.000000000000, 1.000000000000],
        "components": ["A", "B", "Va"],
        "independent_compositions": 2,
        "mol_formula": "A(1-a-b)B(a)Va(b)",
        "origin": [1.000000000000, 0.000000000000, 0.000000000000],
        "param_formula": "a(0.333333-0.333333A+0.666667B-0.333333Va)b(0.333333-0.333333A-0.333333B+0.666667Va)"
    }
    return CompositionConverter.from_dict(data)
