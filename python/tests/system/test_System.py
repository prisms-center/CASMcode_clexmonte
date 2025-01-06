import io
from contextlib import redirect_stdout

import libcasm.clusterography as casmclust
import libcasm.configuration as casmconfig
import libcasm.occ_events as occ_events
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


def test_System_from_dict_1_verbose(
    FCCBinaryVacancy_system_data, session_shared_datadir
):
    f = io.StringIO()
    with redirect_stdout(f):
        system = System.from_dict(
            data=FCCBinaryVacancy_system_data,
            search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
            verbose=True,
        )
    out = f.getvalue()
    assert 'Parsing required "prim"...' in out

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


def test_kmc_System_from_dict_1(
    FCCBinaryVacancy_kmc_system_data, session_shared_datadir
):
    system = System.from_dict(
        data=FCCBinaryVacancy_kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
        verbose=False,
    )
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)

    ## Check basis set cluster info
    assert system.basis_set_keys == ["default"]

    def _check(
        orbits,
        function_to_orbit_index,
        expected_orbit_sizes,
        expected_function_to_orbit_index,
    ):
        assert isinstance(orbits, list)
        assert len(orbits) == len(expected_orbit_sizes)
        orbit_sizes = [len(orbit) for orbit in orbits]
        assert orbit_sizes == expected_orbit_sizes
        for orbit in orbits:
            assert isinstance(orbit, list)
            for cluster in orbit:
                assert isinstance(cluster, casmclust.Cluster)
        assert function_to_orbit_index == expected_function_to_orbit_index

    _check(
        *system.basis_set_cluster_info("default"),
        expected_orbit_sizes=[1, 1, 6, 3],
        expected_function_to_orbit_index=[0, 1, 1, 2, 2, 2, 3, 3, 3],
    )

    ## Check local basis set cluster info
    assert system.local_basis_set_keys == ["A_Va_1NN", "B_Va_1NN"]

    def _check(
        equiv_orbits,
        function_to_orbit_index,
        expected_n_equiv,
        expected_orbit_sizes,
        expected_function_to_orbit_index,
    ):
        assert isinstance(equiv_orbits, list)
        assert len(equiv_orbits) == expected_n_equiv
        for orbits in equiv_orbits:
            assert len(orbits) == len(expected_orbit_sizes)
            orbit_sizes = [len(orbit) for orbit in orbits]
            assert orbit_sizes == expected_orbit_sizes
            for orbit in orbits:
                assert isinstance(orbit, list)
                for cluster in orbit:
                    assert isinstance(cluster, casmclust.Cluster)
        assert function_to_orbit_index == expected_function_to_orbit_index

    _check(
        *system.local_basis_set_cluster_info("A_Va_1NN"),
        expected_n_equiv=6,
        expected_orbit_sizes=[1, 2, 4, 4, 8],
        expected_function_to_orbit_index=[0, 1, 1, 2, 2, 3, 3, 4, 4],
    )
    _check(
        *system.local_basis_set_cluster_info("B_Va_1NN"),
        expected_n_equiv=6,
        expected_orbit_sizes=[1, 2, 4, 4, 8],
        expected_function_to_orbit_index=[0, 1, 1, 2, 2, 3, 3, 4, 4],
    )

    assert system.event_type_names == ["A_Va_1NN", "B_Va_1NN"]

    events = system.events("A_Va_1NN")
    assert len(events) == 6
    for event in events:
        assert isinstance(event, occ_events.OccEvent)

    events = system.events("B_Va_1NN")
    assert len(events) == 6
    for event in events:
        assert isinstance(event, occ_events.OccEvent)


def test_kmc_System_1(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System
    assert isinstance(system, System)
    assert isinstance(system.xtal_prim, xtal.Prim)
    assert isinstance(system.prim, casmconfig.Prim)
    assert isinstance(system.n_dimensions, int)
    assert isinstance(system.composition_converter, CompositionConverter)
    assert isinstance(system.composition_calculator, CompositionCalculator)
    assert isinstance(system.prim_neighbor_list, PrimNeighborList)
