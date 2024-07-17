import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.configuration as casmconfig


def test_make_canonical_initial_state_1(Clex_ZrO_Occ_System):
    """Test default motif"""
    system = Clex_ZrO_Occ_System

    # The mol composition element meaning is determined by the
    # order of components in the composition calculator
    assert system.composition_calculator.components() == ["Zr", "Va", "O"]

    # The param_composition meaning is determined by the origin and end member
    # mol compositions
    assert np.allclose(system.composition_converter.origin(), [2.0, 2.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(0), [2.0, 0.0, 2.0])

    supercells = casmconfig.SupercellSet(prim=system.prim)
    motif_in = casmconfig.Configuration.from_dict(
        data={
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        supercells=supercells,
    )

    calculator = clexmonte.MonteCalculator(
        method="canonical",
        system=system,
    )

    modifying_functions = calculator.modifying_functions
    assert isinstance(modifying_functions, clexmonte.StateModifyingFunctionMap)

    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={"temperature": 300.0},  # no composition -> use motif composition
        dirs="abc",
        min_volume=1000,
        motif=motif_in,
    )

    assert isinstance(initial_state, clexmonte.MonteCarloState)
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super.tolist()
        == [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"],
        np.array([2.0, 1.6666666666666667, 0.3333333333333333]),
    )
    assert motif == motif_in


def test_make_canonical_initial_state_2(Clex_ZrO_Occ_System):
    """Test default motif"""
    system = Clex_ZrO_Occ_System

    # The mol composition element meaning is determined by the
    # order of components in the composition calculator
    assert system.composition_calculator.components() == ["Zr", "Va", "O"]

    # The param_composition meaning is determined by the origin and end member
    # mol compositions
    assert np.allclose(system.composition_converter.origin(), [2.0, 2.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(0), [2.0, 0.0, 2.0])

    supercells = casmconfig.SupercellSet(prim=system.prim)
    motif_in = casmconfig.Configuration.from_dict(
        data={
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        supercells=supercells,
    )

    calculator = clexmonte.MonteCalculator(
        method="canonical",
        system=system,
    )

    modifying_functions = calculator.modifying_functions
    assert isinstance(modifying_functions, clexmonte.StateModifyingFunctionMap)

    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={"temperature": 300.0, "mol_composition": [2.0, 1.7, 0.3]},
        dirs="abc",
        min_volume=1000,
        motif=motif_in,
    )

    assert isinstance(initial_state, clexmonte.MonteCarloState)
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super.tolist()
        == [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"],
        np.array([2.0, 1.6997084548104957, 0.30029154518950435]),
    )
    assert motif == motif_in
