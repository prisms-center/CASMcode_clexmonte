import math

import numpy as np

import libcasm.clexmonte as clexmonte


def test_MonteEventData_1(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    state = clexmonte.MonteCarloState(
        configuration=system.make_default_configuration(
            transformation_matrix_to_super=np.eye(3, dtype="int") * 2,
        ),
        conditions={
            "temperature": 300.0,
            "param_composition": [0.0, 0.0],  # <-one of param/mol composition is needed
        },
    )
    event_data = clexmonte.MonteEventData(calculator=calculator, state=state)

    assert len(event_data.prim_event_list) == 24
    assert len(event_data.event_list) == 24 * 8
    assert isinstance(event_data.event_list.total_rate(), float)
    assert math.isclose(event_data.event_list.total_rate(), 0.0, abs_tol=1e-10)
