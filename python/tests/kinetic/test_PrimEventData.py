import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.xtal as xtal


def test_PrimEventData_1(FCCBinaryVacancy_kmc_System):
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

    # 6 A-Va exchange events, 6 B-Va exchange events, x2 for forward and reverse
    assert len(event_data.prim_event_list) == 24

    # -- check PrimEventData 0 --
    prim_event_data = event_data.prim_event_list[0]
    assert isinstance(prim_event_data, clexmonte.PrimEventData)
    assert prim_event_data.equivalent_index == 0
    assert prim_event_data.event_type_name == "A_Va_1NN"
    assert prim_event_data.event == system.events("A_Va_1NN")[0]
    assert prim_event_data.is_forward is True
    assert prim_event_data.prim_event_index == 0
    assert list(prim_event_data.occ_init) == [0, 2]
    assert list(prim_event_data.occ_final) == [2, 0]
    sites = prim_event_data.sites
    assert len(sites) == 2
    for site in sites:
        assert isinstance(site, xtal.IntegralSiteCoordinate)

    # -- check PrimEventData 1 --
    prim_event_data = event_data.prim_event_list[1]
    assert prim_event_data.equivalent_index == 0
    assert prim_event_data.event == system.events("A_Va_1NN")[0].copy_reverse()
    assert prim_event_data.is_forward is False
    assert prim_event_data.prim_event_index == 1
    assert list(prim_event_data.occ_init) == [2, 0]
    assert list(prim_event_data.occ_final) == [0, 2]

    # -- check PrimEventData 2 --
    prim_event_data = event_data.prim_event_list[2]
    assert isinstance(prim_event_data, clexmonte.PrimEventData)
    assert prim_event_data.equivalent_index == 1
    assert prim_event_data.event_type_name == "A_Va_1NN"
    assert prim_event_data.event == system.events("A_Va_1NN")[1]
    assert prim_event_data.is_forward is True
    assert prim_event_data.prim_event_index == 2
    assert list(prim_event_data.occ_init) == [0, 2]
    assert list(prim_event_data.occ_final) == [2, 0]
    sites = prim_event_data.sites
    assert len(sites) == 2
    for site in sites:
        assert isinstance(site, xtal.IntegralSiteCoordinate)

    # -- check PrimEventData 12 --
    prim_event_data = event_data.prim_event_list[12]
    assert isinstance(prim_event_data, clexmonte.PrimEventData)
    assert prim_event_data.equivalent_index == 0
    assert prim_event_data.event_type_name == "B_Va_1NN"
    assert prim_event_data.event == system.events("B_Va_1NN")[0]
    assert prim_event_data.is_forward is True
    assert prim_event_data.prim_event_index == 12
    assert list(prim_event_data.occ_init) == [1, 2]
    assert list(prim_event_data.occ_final) == [2, 1]
    sites = prim_event_data.sites
    assert len(sites) == 2
    for site in sites:
        assert isinstance(site, xtal.IntegralSiteCoordinate)


def test_PrimEventData_repr(FCCBinaryVacancy_kmc_System):
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
    event_system = system.event_system

    # 6 A-Va exchange events, 6 B-Va exchange events, x2 for forward and reverse
    assert len(event_data.prim_event_list) == 24

    # Test PrimEventData repr
    for i, prim_event_data in enumerate(event_data.prim_event_list):
        assert isinstance(prim_event_data, clexmonte.PrimEventData)

        # Test __repr__
        import io
        from contextlib import redirect_stdout

        f = io.StringIO()
        with redirect_stdout(f):
            print(prim_event_data)
        out = f.getvalue()
        assert "equivalent_index" in out
        assert "event" in out
        assert "event_type_name" in out
        assert "is_forward" in out
        assert "occ_init" in out
        assert "occ_final" in out
        assert "prim_event_index" in out
        assert "sites" in out

        # print(f"--- Prim event: {i} ---")
        # print(prim_event_data)
        # print()

        # Test to_dict
        data = prim_event_data.to_dict(event_system=event_system)
        assert "equivalent_index" in data
        assert "event" in data
        assert "event_type_name" in data
        assert "is_forward" in data
        assert "occ_init" in data
        assert "occ_final" in data
        assert "prim_event_index" in data
        assert "sites" in data

        # print(data)
        # print()
