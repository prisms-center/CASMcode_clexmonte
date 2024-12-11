import libcasm.clexmonte as clexmonte


def test_MonteEventData_1(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    vol = 10**3
    Va_comp = 1 / vol
    state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 400.0,
            # one of param/mol composition is needed
            "mol_composition": [0.8 - Va_comp, 0.200, Va_comp],
        },
        min_volume=vol,
    )

    event_data = clexmonte.MonteEventData(calculator=calculator, state=state)
    assert isinstance(event_data, clexmonte.MonteEventData)

    assert len(event_data.prim_event_list) == 24
    assert len(event_data.event_list) == 12
    assert isinstance(event_data.event_list.total_rate(), float)

    event_data_summary = calculator.event_data_summary()
    assert isinstance(event_data_summary, clexmonte.MonteEventDataSummary)

    data = event_data_summary.to_dict()
    assert isinstance(data, dict)

    summary_str = str(event_data_summary)
    assert isinstance(summary_str, str)
    assert "- Number of events (total) = 12" in summary_str
