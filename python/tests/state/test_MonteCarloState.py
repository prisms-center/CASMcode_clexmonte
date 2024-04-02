import copy

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.monte as monte
from libcasm.clexmonte import (
    MonteCarloState,
)


def test_MonteCarloState_constructor_1(
    FCCBinaryVacancy_prim_config,
):
    mc_state = MonteCarloState(
        configuration=FCCBinaryVacancy_prim_config,
    )
    assert isinstance(mc_state, MonteCarloState)
    assert isinstance(mc_state.configuration, casmconfig.Configuration)
    assert isinstance(mc_state.conditions, monte.ValueMap)
    assert isinstance(mc_state.properties, monte.ValueMap)


def test_MonteCarloState_constructor_2(
    FCCBinaryVacancy_prim_config,
):
    mc_state = MonteCarloState(
        configuration=FCCBinaryVacancy_prim_config,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
    )
    assert isinstance(mc_state, MonteCarloState)
    assert isinstance(mc_state.configuration, casmconfig.Configuration)
    assert isinstance(mc_state.conditions, monte.ValueMap)
    assert isinstance(mc_state.properties, monte.ValueMap)


def test_MonteCarloState_copy_1(
    FCCBinaryVacancy_prim_config,
):
    mc_state = MonteCarloState(
        configuration=FCCBinaryVacancy_prim_config,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
    )
    assert isinstance(mc_state, MonteCarloState)
    assert isinstance(mc_state.configuration, casmconfig.Configuration)
    assert isinstance(mc_state.conditions, monte.ValueMap)
    assert isinstance(mc_state.properties, monte.ValueMap)

    mc_state_2 = copy.copy(mc_state)
    assert mc_state is not mc_state_2
    assert mc_state.configuration == mc_state_2.configuration
    assert not mc_state.conditions.is_mismatched(mc_state_2.conditions)
    assert np.allclose(
        mc_state.conditions.scalar_values["temperature"],
        mc_state_2.conditions.scalar_values["temperature"],
    )
    assert np.allclose(
        mc_state.conditions.vector_values["param_chem_pot"],
        mc_state_2.conditions.vector_values["param_chem_pot"],
    )


def test_MonteCarloState_to_from_dict(
    FCCBinaryVacancy_prim_config,
):
    supercells = casmconfig.SupercellSet(
        prim=FCCBinaryVacancy_prim_config.supercell.prim
    )
    supercell = supercells.add_by_transformation_matrix_to_super(
        np.eye(3, dtype="int") * 2,
    ).supercell
    mc_state = MonteCarloState(
        configuration=casmconfig.Configuration(supercell=supercell),
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
    )
    mc_state.configuration.set_occ(0, 1)

    # to_dict
    data = mc_state.to_dict()
    assert isinstance(data, dict)
    assert "configuration" in data
    assert "dof" in data["configuration"]
    assert "occ" in data["configuration"]["dof"]
    assert data["configuration"]["dof"]["occ"] == [1] + 7 * [0]
    assert "transformation_matrix_to_supercell" in data["configuration"]
    assert data["configuration"]["transformation_matrix_to_supercell"] == [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
    ]
    assert "conditions" in data
    assert np.allclose(data["conditions"]["temperature"], 300.0)
    assert np.allclose(data["conditions"]["param_chem_pot"], [-1.0])

    # from_dict
    mc_state_2 = MonteCarloState.from_dict(
        data=data,
        supercells=supercells,
    )
    assert isinstance(mc_state_2, MonteCarloState)
    assert mc_state.configuration == mc_state_2.configuration
    assert np.allclose(
        mc_state.conditions.scalar_values["temperature"],
        mc_state_2.conditions.scalar_values["temperature"],
    )
    assert np.allclose(
        mc_state.conditions.vector_values["param_chem_pot"],
        mc_state_2.conditions.vector_values["param_chem_pot"],
    )
