import libcasm.configuration as casmconfig
import libcasm.xtal as xtal
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


def test_System_from_dict_1(FCCBinaryVacancy_system_data, session_shared_datadir):
    system = System.from_dict(
        data=FCCBinaryVacancy_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_System_1(FCCBinaryVacancy_System):
    system = FCCBinaryVacancy_System
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_System_from_dict_2(FCCBinaryVacancy_system_data, session_shared_datadir):
    system = System.from_dict(
        data=FCCBinaryVacancy_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_kmc_System_from_dict_1(
    FCCBinaryVacancy_kmc_system_data, session_shared_datadir
):
    system = System.from_dict(
        data=FCCBinaryVacancy_kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_kmc_System_from_dict_2(
    FCCBinaryVacancy_kmc_system_data, session_shared_datadir
):
    system = System.from_dict(
        data=FCCBinaryVacancy_kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)


def test_kmc_System_1(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)
