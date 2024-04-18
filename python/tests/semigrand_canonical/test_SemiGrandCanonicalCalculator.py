import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
import libcasm.monte as monte
import libcasm.xtal as xtal


def test_constructors_1(Clex_ZrO_Occ_System):
    system = Clex_ZrO_Occ_System

    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )
    assert isinstance(calculator, clexmonte.MonteCalculator)
    with pytest.raises(Exception):
        assert isinstance(calculator.potential, clexmonte.MontePotential)
    with pytest.raises(Exception):
        assert isinstance(calculator.state_data, clexmonte.StateData)

    state = clexmonte.MonteCarloState(
        configuration=system.make_default_configuration(
            transformation_matrix_to_super=np.eye(3, dtype="int") * 2,
        ),
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [0.0],
        },
    )
    calculator.set_state_and_potential(state=state)
    assert isinstance(calculator.potential, clexmonte.MontePotential)
    assert isinstance(calculator.state_data, clexmonte.StateData)

    state_data = clexmonte.StateData(
        system=system,
        state=state,
        occ_location=None,
    )
    assert isinstance(state_data, clexmonte.StateData)

    potential = clexmonte.MontePotential(calculator=calculator, state=state)
    assert isinstance(potential, clexmonte.MontePotential)


def test_run_fixture_1(Clex_ZrO_Occ_System, tmp_path):
    """A single run, using a fixture"""
    system = Clex_ZrO_Occ_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )

    # construct default sampling fixture parameters
    thermo = calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )
    print(xtal.pretty_json(thermo.to_dict()))

    # construct the initial state (default configuration)
    initial_state, motif, motif_id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        min_volume=1000,
    )

    # Run
    sampling_fixture = calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=thermo,
    )
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    pytest.helpers.validate_summary_file(summary_file=summary_file, expected_size=1)


def test_run_fixture_2(Clex_ZrO_Occ_System, tmp_path):
    """A single run, using a fixture, for a range of param_chem_pot"""
    system = Clex_ZrO_Occ_System

    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )

    # construct default sampling fixture parameters
    thermo = calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )

    # set lower convergence level for potential_energy
    thermo.converge(quantity="potential_energy", abs=2e-3)

    # set lower convergence level for param_composition("a")
    thermo.converge(quantity="param_composition", abs=2e-3, component_name=["a"])

    # construct the initial state (default configuration)
    state, motif, motif_id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        min_volume=1000,
    )

    # Run several, w/ dependent runs
    x_list = np.arange(-4.0, 0.01, step=0.5)
    for x in x_list:
        state.conditions.vector_values["param_chem_pot"] = [x]
        sampling_fixture = calculator.run_fixture(
            state=state,
            sampling_fixture_params=thermo,
        )
        assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    pytest.helpers.validate_summary_file(
        summary_file=summary_file, expected_size=len(x_list)
    )


def test_run_1(Clex_ZrO_Occ_System, tmp_path):
    """A single run, using RunManager"""
    system = Clex_ZrO_Occ_System

    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="semigrand_canonical",
        system=system,
    )

    # construct default sampling fixture parameters
    thermo = calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )

    # construct RunManager
    run_manager = clexmonte.RunManager(
        engine=monte.RandomNumberEngine(),
        sampling_fixture_params=[thermo],
        global_cutoff=True,
    )

    # construct the initial state (default configuration)
    state, motif, motif_id = clexmonte.make_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        min_volume=1000,
    )

    # Run several, w/ dependent runs
    x_list = np.arange(-4.0, 0.01, step=0.5)
    for x in x_list:
        state.conditions.vector_values["param_chem_pot"] = [x]
        run_manager = calculator.run(
            state=state,
            run_manager=run_manager,
        )
        assert isinstance(run_manager, clexmonte.RunManager)
    assert "thermo" in run_manager.sampling_fixture_labels
    sampling_fixture = run_manager.sampling_fixture("thermo")
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    pytest.helpers.validate_summary_file(
        summary_file=summary_file, expected_size=len(x_list)
    )
