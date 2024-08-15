import copy

import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.configuration as casmconfig
import libcasm.enumerate as casmenum


def test_make_initial_state_1(Clex_ZrO_Occ_System):
    """Test default motif"""
    system = Clex_ZrO_Occ_System

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
        motifs=None,
    )
    assert isinstance(initial_state, clexmonte.MonteCarloState)
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super
        == np.eye(3, dtype="int") * 10
    ).all()
    assert (
        motif.supercell.transformation_matrix_to_super == np.eye(3, dtype="int")
    ).all()
    assert (motif.occupation == np.array([0, 0, 0, 0], dtype="int")).all()
    assert id == "default"


def test_make_initial_state_2a(Clex_ZrO_Occ_System):
    """Test motif"""
    system = Clex_ZrO_Occ_System

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
        method="semigrand_canonical",
        system=system,
    )
    initial_state, motif_out, id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        dirs="abc",
        min_volume=1000,
        motifs=[motif_in],
    )
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super.tolist()
        == [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ]
    )
    assert motif_out == motif_in
    assert id == 0


def test_make_initial_state_2b(Clex_ZrO_Occ_System):
    """Test motif w/ transformation_matrix_to_super"""
    system = Clex_ZrO_Occ_System

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
        method="semigrand_canonical",
        system=system,
    )
    initial_state, motif_out, id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        dirs="abc",
        min_volume=1000,
        transformation_matrix_to_super=np.array(
            [
                [2, 1, 2],
                [1, 2, 2],
                [0, 0, 2],
            ],
            dtype="int",
        ),
        motifs=[motif_in],
    )
    assert (
        initial_state.configuration.supercell.transformation_matrix_to_super.tolist()
        == [
            [12, 6, 12],
            [6, 12, 12],
            [0, 0, 12],
        ]
    )
    assert motif_out == motif_in
    assert id == 0


def test_make_initial_state_3(Clex_ZrO_Occ_System):
    """Test motifs=ConfigurationSet"""
    system = Clex_ZrO_Occ_System

    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )

    ###
    supercells = casmconfig.SupercellSet(prim=system.prim)
    configurations = casmconfig.ConfigurationSet()
    enum = casmenum.ConfigEnumAllOccupations(prim=system.prim, supercell_set=supercells)
    for config in enum.by_supercell(
        max=6,
        skip_non_canonical=True,
        skip_non_primitive=True,
    ):
        configurations.add(config)
    ###

    initial_state_list = []
    motif_list = []
    id_list = []

    x_range = np.arange(-4.0, 0.01, step=0.1)
    for x in x_range:
        initial_state, motif, id = clexmonte.make_initial_state(
            calculator=calculator,
            conditions={
                "temperature": 300.0,
                "param_chem_pot": [x],
            },
            dirs="abc",
            min_volume=1000,
            motifs=configurations,
        )
        if id not in id_list:
            initial_state_list.append(initial_state)
            motif_list.append(motif)
            id_list.append(id)

    T_data = [
        x.configuration.supercell.transformation_matrix_to_super.tolist()
        for x in initial_state_list
    ]
    print(T_data)
    assert T_data == [
        [
            [10, 0, 0],
            [0, 10, 0],
            [0, 0, 10],
        ],
        [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ],
        [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ],
        [
            [0, 14, -7],
            [0, 7, 7],
            [7, 0, 0],
        ],
    ]

    motif_data = [x.to_dict() for x in motif_list]
    assert motif_data == [
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0]},
            "supercell_name": "SCEL1_1_1_1_0_0_0",
            "transformation_matrix_to_supercell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1]},
            "supercell_name": "SCEL3_3_1_1_0_0_2",
            "transformation_matrix_to_supercell": [[0, 2, -1], [0, 1, 1], [1, 0, 0]],
        },
    ]
    assert id_list == [
        "SCEL1_1_1_1_0_0_0/0",
        "SCEL3_3_1_1_0_2_2/0",
        "SCEL3_3_1_1_0_2_2/5",
        "SCEL3_3_1_1_0_0_2/8",
    ]


def test_make_initial_state_4(Clex_ZrO_Occ_System):
    """Test motifs=ConfigurationSet"""
    system = Clex_ZrO_Occ_System
    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )

    ###
    supercells = casmconfig.SupercellSet(prim=system.prim)
    configuration_list = []
    enum = casmenum.ConfigEnumAllOccupations(prim=system.prim, supercell_set=supercells)
    for config in enum.by_supercell(
        max=4,
        skip_non_canonical=False,
        skip_non_primitive=False,
    ):
        configuration_list.append(copy.copy(config))
    ###
    print(len(configuration_list))

    initial_state_list = []
    motif_list = []
    id_list = []

    x_range = np.arange(-4.0, 0.01, step=0.1)
    for x in x_range:
        initial_state, motif, id = clexmonte.make_initial_state(
            calculator=calculator,
            conditions={
                "temperature": 300.0,
                "param_chem_pot": [x],
            },
            dirs="abc",
            min_volume=1000,
            motifs=configuration_list,
        )
        if id not in id_list:
            initial_state_list.append(initial_state)
            motif_list.append(motif)
            id_list.append(id)

    T_data = [
        x.configuration.supercell.transformation_matrix_to_super.tolist()
        for x in initial_state_list
    ]
    assert T_data == [
        [
            [10, 0, 0],
            [0, 10, 0],
            [0, 0, 10],
        ],
        [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ],
        [
            [14, 7, 7],
            [7, 14, 7],
            [0, 0, 7],
        ],
        [
            [0, 14, -7],
            [0, 7, 7],
            [7, 0, 0],
        ],
    ]

    motif_data = [x.to_dict() for x in motif_list]
    assert motif_data == [
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0]},
            "supercell_name": "SCEL1_1_1_1_0_0_0",
            "transformation_matrix_to_supercell": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_2_2",
            "transformation_matrix_to_supercell": [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
        },
        {
            "basis": "standard",
            "dof": {"occ": [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0]},
            "supercell_name": "SCEL3_3_1_1_0_0_2",
            "transformation_matrix_to_supercell": [[0, 2, -1], [0, 1, 1], [1, 0, 0]],
        },
    ]
    assert id_list == [0, 245, 256, 130]
