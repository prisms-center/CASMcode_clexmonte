import contextlib

import numpy as np
import pytest

import libcasm.clexmonte as clexmonte

CalculatorTestRunner = pytest.helpers.CalculatorTestRunner


def setup_1(runner: CalculatorTestRunner):
    """Adds the following:

    - A custom event state calculation function for the A_Va_1NN and B_Va_1NN events,
      which in turn just prints a message, calls the default event state calculation
      function, and increments a counter.
    - A custom selected event function, which just prints a message and increments a
      counter.

    runner.kinetics_params: SamplingFixtureParams
        Sampling fixture parameters:
        - sample_by_step()
        - set_max_count(100)

    runner.n_event_calculations: int
        The count of the number of calls of the custom event state calculation.

    runner.n_selected_events: int
        The count of the number of calls of the custom selected event function.

    Parameters
    ----------
    runner: CalculatorTestRunner
        The runner object with the MonteCalculator.

    Returns
    -------
    runner: CalculatorTestRunner
        The runner object with the added attributes.
    """

    kinetics_params = runner.calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(runner.output_dir / "kinetics"),
        write_observations=False,
    )
    kinetics_params.sample_by_pass()
    kinetics_params.set_max_count(100)
    runner.kinetics_params = kinetics_params

    calculator = runner.calculator
    # calculator.collect("selected_event.by_equivalent_index")
    # calculator.collect("dE_activated.by_type", bin_width=0.01)

    # A custom event state calculation function which just increments a counter, prints
    # a message, and calls the default event state calculation
    runner.n_event_calculations = 0
    runner.n_encountered_abnormal = 0

    def event_state_calculation_f(
        state: clexmonte.EventState,  # <-- this is a mutable reference
        event_state_calculator: clexmonte.EventStateCalculator,
    ):
        runner.n_event_calculations += 1

        # print(f"CUSTOM EVENT STATE CALCULATION {runner.n_event_calculations}")
        #
        # prim_event_data = event_state_calculator.curr_prim_event_data
        # print(
        #     "event_type_name:",
        #     prim_event_data.event_type_name,
        #     "unitcell_index:",
        #     event_state_calculator.curr_unitcell_index,
        #     "equivalent_index:",
        #     prim_event_data.equivalent_index,
        # )
        event_state_calculator.set_default_event_state(state)
        if not state.is_normal:
            runner.n_encountered_abnormal += 1
            # print(f"- NOT NORMAL {runner.n_encountered_abnormal}")
        # print(
        #     "dE_activated:",
        #     state.dE_activated,
        #     "dE_final:",
        #     state.dE_final,
        #     "is_normal:",
        #     state.is_normal,
        # )
        # print(
        #     "rate:",
        #     state.rate,
        # )
        # print()
        # sys.stdout.flush()

    calculator.event_data.set_custom_event_state_calculation(
        event_type_name="A_Va_1NN",
        function=event_state_calculation_f,
    )

    calculator.event_data.set_custom_event_state_calculation(
        event_type_name="B_Va_1NN",
        function=event_state_calculation_f,
    )

    # A custom selected event function which just increments a counter and prints
    # a message
    runner.n_selected_events = 0

    def print_selected_event_f():
        runner.n_selected_events += 1
        # print(f"Event was selected {runner.n_selected_events}")

    calculator.add_generic_function(
        name="print_selected_event_f",
        description="Print information about the selected event state",
        requires_event_state=False,
        function=print_selected_event_f,
        order=0,
    )
    calculator.evaluate("print_selected_event_f")

    return runner


def add_initial_state(
    runner: CalculatorTestRunner,
    temperature: float,
    mol_composition: np.ndarray,
    min_volume: int = 1000,
):
    """Adds the following:

    runner.initial_state: MonteCarloState
        The initial state is set to a canonical state with
        temperature=1200.0 and mol_composition=[0.899, 0.1, 0.001],
        and volume=1000.

    runner.state: MonteCarloState
        A copy of `runner.initial_state`.

    runner.motif: Configuration
        The initial state's motif (default configuration)

    Parameters
    ----------
    runner: CalculatorTestRunner
        The runner object with the MonteCalculator.

    Returns
    -------
    runner: CalculatorTestRunner
        The runner object with the added attributes.
    """
    calculator = runner.calculator
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": temperature,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": mol_composition,
        },
        min_volume=min_volume,
    )
    runner.initial_state = initial_state
    runner.state = initial_state.copy()
    runner.motif = motif
    return runner


def run_test(
    runner: CalculatorTestRunner,
    temperature: float,
    mol_composition: list,
    seed: int,
    expect_set_event_data_exception: bool,
    expect_run_fixture_exception: bool,
    expected_n_event_calculations: int,
    expected_n_encountered_abnormal: dict,
    expected_n_selected_abnormal: dict,
    expected_n_write_encountered: int,
    expected_n_write_selected: int,
):
    """Runs the test with the given parameters.

    Parameters
    ----------
    runner: CalculatorTestRunner
        The runner object with the MonteCalculator.

    params: dict
        The parameters to set in the calculator.

    temperature: float
        The temperature to set in the initial state.

    mol_composition: list
        The mol composition to set in the initial state.

    seed: int
        The seed to set in the calculator.

    expect_set_event_data_exception: bool, optional
        Whether to expect an exception when calling `set_event_data`.
        Default is True.

    expect_run_fixture_exception: bool, optional
        Whether to expect an exception when calling `run_fixture`.
        Default is True.

    expected_n_event_calculations: int, optional
        The expected number of event calculations.
        Default is 100.

    expected_n_encountered_abnormal: dict, optional
        The expected number of encountered abnormal events, by event type, from
        `runner.calculator.event_data.n_encountered_abnormal`.
        Default is {}.

    expected_n_selected_abnormal: dict, optional
        The expected number of selected abnormal events, by event type, from
        `runner.calculator.event_data.n_selected_abnormal`.
        Default is {}.
    """
    runner = setup_1(runner)
    runner.calculator.engine.seed(seed)
    runner = add_initial_state(
        runner=runner,
        temperature=temperature,
        mol_composition=mol_composition,
    )
    runner.calculator.set_state_and_potential(runner.state)
    runner.calculator.make_occ_location()

    print(f"Test output dir: {runner.output_dir}\n")

    if expect_set_event_data_exception:
        # Calculate event rates - expect an exception due to an abnormal event
        with pytest.raises(Exception):
            runner.calculator.set_event_data()

    elif expect_run_fixture_exception:
        # Run - expect an exception due to an abnormal event
        with pytest.raises(Exception):
            runner.calculator.run_fixture(
                state=runner.state,
                sampling_fixture_params=runner.kinetics_params,
            )

    else:
        # Run - expect an exception due to an abnormal event
        runner.calculator.run_fixture(
            state=runner.state,
            sampling_fixture_params=runner.kinetics_params,
        )

    event_data = runner.calculator.event_data

    assert runner.n_event_calculations == expected_n_event_calculations
    assert event_data.n_encountered_abnormal == expected_n_encountered_abnormal
    assert event_data.n_selected_abnormal == expected_n_selected_abnormal

    # --- Test the abnormal_events jsonl files were written ---

    if expected_n_write_encountered or expected_n_write_selected:
        assert runner.output_dir.exists()

        for event_type, expected_count in [
            ("encountered_abnormal_events.jsonl", expected_n_write_encountered),
            ("selected_abnormal_events.jsonl", expected_n_write_selected),
        ]:
            file = runner.output_dir / event_type
            if expected_count != 0:
                assert file.exists()

                with open(file, "r") as f:
                    lines = 0
                    for line in f:
                        lines += 1

                assert lines == expected_count

            else:
                assert not file.exists()

    # --- Test the read_abnormal_events function ---

    def _summarize(local_configurations, event_data):
        for event_type_name in local_configurations:
            _event_data = event_data[event_type_name]
            print(f"Event type name: {event_type_name}")
            if len(_event_data) == 0:
                print("- No abnormal events")
            else:
                for i, entry in enumerate(_event_data):
                    s = entry.get("event_state")
                    print(
                        f"- {i}: ",
                        f"Ekra={s['Ekra']:.6f}",
                        f"dE_final={s['dE_final']:.6f}",
                        f"dE_activated= {s['dE_activated']:.6f}",
                    )
            print()

    # which=="all"
    calculator = runner.calculator
    n, local_configurations, event_data = calculator.read_abnormal_events()
    if expected_n_write_encountered + expected_n_write_selected == 0:
        assert n is None
    else:
        assert n == expected_n_write_encountered + expected_n_write_selected
    # print("all:")
    # _summarize(local_configurations, event_data)

    # which=="encountered"
    n, local_configurations, event_data = calculator.read_abnormal_events(
        which="encountered",
    )
    if expected_n_write_encountered == 0:
        assert n is None
    else:
        assert n == expected_n_write_encountered
    # print("encountered:")
    # _summarize(local_configurations, event_data)

    # which=="selected"
    n, local_configurations, event_data = calculator.read_abnormal_events(
        which="selected",
    )
    if expected_n_write_selected == 0:
        assert n is None
    else:
        assert n == expected_n_write_selected
    # print("selected:")
    # _summarize(local_configurations, event_data)


def test_abnormal_event_handling_1a(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test the default handling of abnormal events.

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """
    with contextlib.chdir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    # "abnormal_event_handling": None,
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=699,
            expect_set_event_data_exception=True,
            expect_run_fixture_exception=True,
            expected_n_event_calculations=3,
            expected_n_encountered_abnormal={"B_Va_1NN": 1},
            expected_n_selected_abnormal={},
            expected_n_write_encountered=1,
            expected_n_write_selected=0,
        )


def test_abnormal_event_handling_1b(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test only throw on selected abnormal events

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """

    with contextlib.chdir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                    }
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=699,
            expect_set_event_data_exception=False,
            expect_run_fixture_exception=True,
            expected_n_event_calculations=98423,
            expected_n_encountered_abnormal={"B_Va_1NN": 209},
            expected_n_selected_abnormal={"B_Va_1NN": 1},
            expected_n_write_encountered=27,
            expected_n_write_selected=1,
        )


def test_abnormal_event_handling_1c(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test only throw on selected abnormal events

    with lower temperature - encountered but did not select an abnormal event


    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """
    with contextlib.chdir(tmp_path):
        # Test only throw on selected abnormal events:
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "debug",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                    },
                },
                output_dir=tmp_path / "output",
            ),
            temperature=300.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=699,
            expect_set_event_data_exception=False,
            expect_run_fixture_exception=False,
            expected_n_event_calculations=1300000,
            expected_n_encountered_abnormal={"B_Va_1NN": 1},
            expected_n_selected_abnormal={},
            expected_n_write_encountered=1,
            expected_n_write_selected=0,
        )


def test_abnormal_event_handling_1d(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test no throw on abnormal events

    Seed = 699, mol_composition=[0.799, 0.2, 0.001],
    for an initial state that has an abnormal event.
    """

    with contextlib.chdir(tmp_path):
        run_test(
            runner=CalculatorTestRunner(
                system=FCCBinaryVacancy_kmc_System_2,
                method="kinetic",
                params={
                    # "verbosity": "quiet",
                    "abnormal_event_handling": {
                        "encountered_events": {
                            "throw": False,
                        },
                        "selected_events": {
                            "throw": False,
                        },
                    }
                },
                output_dir=tmp_path / "output",
            ),
            temperature=1200.0,
            mol_composition=[0.799, 0.2, 0.001],
            seed=699,
            expect_set_event_data_exception=False,
            expect_run_fixture_exception=False,
            expected_n_event_calculations=1300000,
            expected_n_encountered_abnormal={"B_Va_1NN": 39981},
            expected_n_selected_abnormal={"B_Va_1NN": 48},
            expected_n_write_encountered=100,
            expected_n_write_selected=10,
        )
