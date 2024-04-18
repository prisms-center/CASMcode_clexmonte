import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.configuration as casmconfig


def mol_composition(
    system: clexmonte.System,
    config: casmconfig.Configuration,
):
    return system.composition_calculator.mean_num_each_component(
        config.dof_values.occupation()
    )


def test_enforce_composition_1(Clex_ZrO_Occ_System):
    """Test default motif"""
    system = Clex_ZrO_Occ_System

    supercells = casmconfig.SupercellSet(
        prim=system.prim,
    )
    motif = casmconfig.Configuration.from_dict(
        data={
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        supercells=supercells,
    )

    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )
    initial_state, motif, id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [0.0],
        },
        dirs="abc",
        min_volume=1000,
        motifs=[motif],
    )

    # n = [n_Zr, n_Va, n_O]
    n_init = mol_composition(system, initial_state.configuration)
    assert np.allclose(n_init, np.array([2.0, 5.0 / 3.0, 1.0 / 3.0]))

    target_mol_composition = np.array([2.0, 1.8, 0.2])
    clexmonte.enforce_composition(
        state=initial_state,
        target_mol_composition=target_mol_composition,
        system=system,
    )

    n_final = mol_composition(system, initial_state.configuration)
    assert np.allclose(n_final, np.array([2.0, 1852 / 1029, 206 / 1029]))
