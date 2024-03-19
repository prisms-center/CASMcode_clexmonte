import copy
import numpy as np

from libcasm.clexmonte import (
    MonteCarloConfiguration,
)
from libcasm.clexulator import (
    ConfigDoFValues,
    make_default_config_dof_values,
)
import libcasm.configuration as casmconfig


def test_MonteCarloConfiguration_constructor_1(
    FCCBinaryVacancy_xtal_prim,
):
    T = np.eye(3, dtype="int")
    n_unitcells = int(round(np.linalg.det(T)))
    mc_config = MonteCarloConfiguration(
        transformation_matrix_to_super=T,
        dof_values=make_default_config_dof_values(
            xtal_prim=FCCBinaryVacancy_xtal_prim,
            n_unitcells=n_unitcells,
        ),
    )
    assert isinstance(mc_config, MonteCarloConfiguration)
    assert isinstance(mc_config.transformation_matrix_to_super, np.ndarray)
    assert isinstance(mc_config.dof_values, ConfigDoFValues)


def test_MonteCarloConfiguration_copy_1(
    FCCBinaryVacancy_xtal_prim,
):
    T = np.eye(3, dtype="int")
    n_unitcells = int(round(np.linalg.det(T)))
    mc_config = MonteCarloConfiguration(
        transformation_matrix_to_super=T,
        dof_values=make_default_config_dof_values(
            xtal_prim=FCCBinaryVacancy_xtal_prim,
            n_unitcells=n_unitcells,
        ),
    )

    x = mc_config.dof_values
    assert x is mc_config.dof_values

    y = copy.copy(mc_config.dof_values)
    assert y is not mc_config.dof_values

    mc_config_2 = mc_config
    assert mc_config_2 is mc_config
    assert x is mc_config_2.dof_values
    assert y is not mc_config_2.dof_values

    mc_config_3 = copy.copy(mc_config)
    assert mc_config_3 is not mc_config
    assert x is not mc_config_3.dof_values
    assert y is not mc_config_3.dof_values


def test_MonteCarloConfiguration_to_from_config_with_sym_info_1(
    FCCBinaryVacancy_xtal_prim,
):
    T = np.eye(3, dtype="int") * 2

    prim = casmconfig.Prim(FCCBinaryVacancy_xtal_prim)
    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=T,
    )
    config_with_sym_info = casmconfig.Configuration(supercell)
    config_with_sym_info.set_occ(0, 1)

    # from_config_with_sym_info
    mc_config = MonteCarloConfiguration.from_config_with_sym_info(config_with_sym_info)

    assert isinstance(mc_config, MonteCarloConfiguration)
    assert np.allclose(mc_config.transformation_matrix_to_super, T)
    assert np.allclose(
        mc_config.dof_values.occupation(), config_with_sym_info.occupation
    )

    # to_config_with_sym_info (prim only)
    config_with_sym_info_2 = mc_config.to_config_with_sym_info(prim=prim)
    assert isinstance(config_with_sym_info_2, casmconfig.Configuration)
    assert np.allclose(
        config_with_sym_info_2.supercell.transformation_matrix_to_super, T
    )
    assert np.allclose(
        config_with_sym_info_2.occupation, config_with_sym_info.occupation
    )

    # to_config_with_sym_info (using SupercellSet)
    supercells = casmconfig.SupercellSet(prim)
    assert len(supercells) == 0

    config_with_sym_info_3 = mc_config.to_config_with_sym_info(supercells=supercells)
    assert len(supercells) == 1
    assert isinstance(config_with_sym_info_3, casmconfig.Configuration)
    assert np.allclose(
        config_with_sym_info_3.supercell.transformation_matrix_to_super, T
    )
    assert np.allclose(
        config_with_sym_info_3.occupation, config_with_sym_info.occupation
    )


def test_MonteCarloConfiguration_to_from_dict(
    FCCBinaryVacancy_xtal_prim,
):
    T = np.eye(3, dtype="int") * 2
    n_unitcells = int(round(np.linalg.det(T)))
    mc_config = MonteCarloConfiguration(
        transformation_matrix_to_super=T,
        dof_values=make_default_config_dof_values(
            xtal_prim=FCCBinaryVacancy_xtal_prim,
            n_unitcells=n_unitcells,
        ),
    )
    mc_config.dof_values.set_occ(0, 1)

    # to_dict
    data = mc_config.to_dict()
    assert isinstance(data, dict)
    assert "dof" in data
    assert "occ" in data["dof"]
    assert data["dof"]["occ"] == [1] + 7 * [0]
    assert "transformation_matrix_to_supercell" in data
    assert data["transformation_matrix_to_supercell"] == [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
    ]

    # from_dict
    mc_config_2 = MonteCarloConfiguration.from_dict(data)
    assert isinstance(mc_config_2, MonteCarloConfiguration)
    assert np.allclose(mc_config_2.transformation_matrix_to_super, T)
    assert np.allclose(
        mc_config_2.dof_values.occupation(), mc_config.dof_values.occupation()
    )
