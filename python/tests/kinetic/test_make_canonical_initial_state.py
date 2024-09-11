import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.configuration as casmconfig


def test_make_canonical_initial_state_1(FCCBinaryVacancy_kmc_System):
    """Test default motif"""
    system = FCCBinaryVacancy_kmc_System

    # The mol composition element meaning is determined by the
    # order of components in the composition calculator
    assert system.composition_calculator.components() == ["A", "B", "Va"]

    # The param_composition meaning is determined by the origin and end member
    # mol compositions
    assert np.allclose(system.composition_converter.origin(), [1.0, 0.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(0), [0.0, 1.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(1), [0.0, 0.0, 1.0])

    supercells = casmconfig.SupercellSet(prim=system.prim)
    motif_in = casmconfig.Configuration.from_dict(
        data={
            "dof": {"occ": [0]},
            "supercell_name": "SCEL1_1_1_1_0_0_0",
            "transformation_matrix_to_supercell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        },
        supercells=supercells,
    )

    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    modifying_functions = calculator.modifying_functions
    assert isinstance(modifying_functions, clexmonte.StateModifyingFunctionMap)

    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "mol_composition": [899 / 1000, 100 / 1000, 1 / 1000],
        },
        dirs="abc",
        min_volume=1000,
        motif=motif_in,
    )

    assert isinstance(initial_state, clexmonte.MonteCarloState)
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super.tolist()
        == [
            [10, 0, 0],
            [0, 10, 0],
            [0, 0, 10],
        ]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"],
        np.array([0.899, 0.1, 0.001]),
    )
    assert motif == motif_in
