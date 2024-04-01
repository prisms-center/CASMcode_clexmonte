import numpy as np
import pytest

import libcasm.configuration as casmconfig


@pytest.fixture
def FCCBinaryVacancy_prim_config(FCCBinaryVacancy_xtal_prim):
    prim = casmconfig.Prim(FCCBinaryVacancy_xtal_prim)

    supercell = casmconfig.Supercell(
        prim=prim,
        transformation_matrix_to_super=np.eye(3, dtype="int"),
    )
    return casmconfig.Configuration(
        supercell=supercell,
    )
