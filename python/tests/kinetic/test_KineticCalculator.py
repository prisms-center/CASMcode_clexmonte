import io
from contextlib import redirect_stdout

import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
import libcasm.monte.events as monte_events
import libcasm.monte.sampling as sampling
import libcasm.xtal as xtal
from libcasm.local_configuration import (
    LocalConfiguration,
)


def test_constructors_1(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System

    # The mol composition element meaning is determined by the
    # order of components in the composition calculator
    assert system.composition_calculator.components() == ["A", "B", "Va"]

    # The param_composition meaning is determined by the origin and end member
    # mol compositions
    assert np.allclose(system.composition_converter.origin(), [1.0, 0.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(0), [0.0, 1.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(1), [0.0, 0.0, 1.0])

    print("-- Construct MonteCalculator --")
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )
    assert isinstance(calculator, clexmonte.MonteCalculator)
    with pytest.raises(Exception):
        assert isinstance(calculator.potential, clexmonte.MontePotential)
    with pytest.raises(Exception):
        assert isinstance(calculator.state_data, clexmonte.StateData)
    assert isinstance(calculator.event_data, clexmonte.MonteEventData)

    # default configuration is occupied by A: [1.0, 0.0, 0.0], which corresponds
    # to the origin composition as defined in the system's composition axes
    print("-- Construct MonteCarloState --")
    state = clexmonte.MonteCarloState(
        configuration=system.make_default_configuration(
            transformation_matrix_to_super=np.eye(3, dtype="int") * 2,
        ),
        conditions={
            "temperature": 300.0,
            "param_composition": [
                0.0,
                0.0,
            ],  # <-one of param/mol composition is needed
            # "mol_composition": [1.0, 0.0, 0.0],
        },
    )

    print("-- Composition checks --")
    composition_calculator = system.composition_calculator
    composition_converter = system.composition_converter
    mol_composition = composition_calculator.mean_num_each_component(
        state.configuration.occupation
    )
    assert np.allclose(mol_composition, [1.0, 0.0, 0.0])
    # assert np.allclose(
    #     mol_composition, state.conditions.vector_values["mol_composition"]
    # )
    param_composition = composition_converter.param_composition(mol_composition)
    assert np.allclose(param_composition, [0.0, 0.0])
    assert np.allclose(
        param_composition, state.conditions.vector_values["param_composition"]
    )

    ## Set state data and potential
    print("-- Set state and potential --")
    calculator.set_state_and_potential(state=state)
    assert isinstance(calculator.potential, clexmonte.MontePotential)
    assert isinstance(calculator.state_data, clexmonte.StateData)
    assert isinstance(calculator.event_data, clexmonte.MonteEventData)

    ## StateData constructor
    print("-- Construct StateData --")
    state_data = clexmonte.StateData(
        system=system,
        state=state,
        occ_location=None,
    )
    assert isinstance(state_data, clexmonte.StateData)

    ## MontePotential constructor
    print("-- Construct MontePotential --")
    potential = clexmonte.MontePotential(calculator=calculator, state=state)
    assert isinstance(potential, clexmonte.MontePotential)

    # Set event data (raise if no occ_location)
    print("-- Check set_event_data w/no occ_location raises --")
    with pytest.raises(Exception):
        # no occ_location set
        calculator.set_event_data()

    ## Make occ_location
    print("-- Make occupant location tracker --")
    occ_location = calculator.make_occ_location()
    assert isinstance(occ_location, monte_events.OccLocation)
    assert isinstance(calculator.event_data, clexmonte.MonteEventData)

    ## Set event data (raise if no allowed events)
    # event_data = clexmonte.MonteEventData(calculator=calculator, state=state)
    print("-- Set event data --")
    calculator.set_event_data()
    event_data = calculator.event_data
    assert isinstance(event_data, clexmonte.MonteEventData)
    assert len(event_data.event_list) == 0

    ## MonteEventData constructor
    print("-- Construct MonteEventData --")
    event_data = clexmonte.MonteEventData(calculator=calculator, state=state)
    assert isinstance(event_data, clexmonte.MonteEventData)
    assert len(event_data.event_list) == 0


def test_event_data_1(FCCBinaryVacancy_kmc_System):
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

    # Test PrimEventData repr
    for i, prim_event_data in enumerate(event_data.prim_event_list):
        assert isinstance(prim_event_data, clexmonte.PrimEventData)
        print(f"--- Prim event: {i} ---")
        print(prim_event_data)
        print()

        f = io.StringIO()
        with redirect_stdout(f):
            print(prim_event_data)
        out = f.getvalue()
        assert "sites" in out


def test_run_fixture_1_sample_state_data(FCCBinaryVacancy_kmc_System, tmp_path):
    """A single run, using a fixture"""
    system = FCCBinaryVacancy_kmc_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )
    calculator.collect("selected_event.by_equivalent_index")

    def print_step_f():
        kinetics_data = calculator.kinetics_data
        assert isinstance(kinetics_data, clexmonte.KineticsData)

        # Check unique_atom_id type
        assert isinstance(kinetics_data.unique_atom_id, clexmonte.LongVector)

        # Check prev_unique_atom_id type
        prev_unique_atom_id = kinetics_data.prev_unique_atom_id
        assert "kinetics" in prev_unique_atom_id
        assert isinstance(prev_unique_atom_id["kinetics"], clexmonte.LongVector)

        # Check atom_name_index_list type
        assert isinstance(kinetics_data.atom_name_index_list, clexmonte.LongVector)

        # Test returning a numpy array with the step number
        fixture = kinetics_data.sampling_fixture
        # print("step:", step)
        return np.array([fixture.n_step])

    print_step_sampling_f = sampling.StateSamplingFunction(
        name="print_step",
        description="test",
        shape=[],
        function=print_step_f,
        component_names=["step"],
    )

    def json_step_f():
        fixture = calculator.kinetics_data.sampling_fixture
        json_step = {"step": fixture.n_step}
        # print("json_step:", json_step)
        return json_step

    json_step_sampling_f = sampling.jsonStateSamplingFunction(
        name="json_step",
        description="test",
        function=json_step_f,
    )

    # construct default sampling fixture parameters
    kinetics = calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(output_dir),
        write_observations=True,
    )
    kinetics.add_sampling_function(print_step_sampling_f)
    kinetics.add_json_sampling_function(json_step_sampling_f)
    kinetics.sample_by_step()
    kinetics.sample("print_step")
    kinetics.sample("json_step")
    kinetics.do_not_sample("clex.formation_energy")
    kinetics.sample("clex.formation_energy")

    kinetics.clear_cutoffs()
    # kinetics.set_min_count(1000)
    kinetics.set_max_count(1000)

    kinetics.converge("clex.formation_energy", abs=1e-3)
    kinetics.do_not_converge("clex.formation_energy")

    print(xtal.pretty_json(kinetics.to_dict()))

    # construct the initial state:
    # start from the default configuration
    # and modify to match mol_composition=[0.899, 0.1, 0.001]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )
    composition_calculator = system.composition_calculator
    assert (
        composition_calculator.num_each_component(
            initial_state.configuration.occupation
        )
        == np.array([899, 100, 1], dtype="int")
    ).all()

    # make_canonical_initial_state should have set the mol_composition
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"], [0.899, 0.1, 0.001]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["param_composition"], [0.1, 0.001]
    )

    # Run
    sampling_fixture = calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=kinetics,
    )
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    # print(sampling_fixture.results.json_samplers["json_step"].values())

    pytest.helpers.validate_summary_file(
        summary_file=summary_file,
        expected_size=1,
        is_canonical=True,
        is_requested_convergence=False,
    )


def test_print_selected_event_functions(FCCBinaryVacancy_kmc_System):
    system = FCCBinaryVacancy_kmc_System

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    calculator.print_selected_event_functions()

    f = io.StringIO()
    with redirect_stdout(f):
        calculator.print_selected_event_functions()
    out = f.getvalue()
    assert "local_orbit_composition.A_Va_1NN-1:" in out
    assert "local_orbit_composition.A_Va_1NN-2:" in out
    assert "local_orbit_composition.A_Va_1NN-3:" in out
    assert "local_orbit_composition.A_Va_1NN-4:" in out
    assert "local_orbit_composition.A_Va_1NN-all:" in out
    assert "local_orbit_composition.A_Va_1NN-all-combined:" in out
    assert "local_orbit_composition.B_Va_1NN-1:" in out
    assert "local_orbit_composition.B_Va_1NN-2:" in out
    assert "local_orbit_composition.B_Va_1NN-3:" in out
    assert "local_orbit_composition.B_Va_1NN-4:" in out
    assert "local_orbit_composition.B_Va_1NN-all:" in out
    assert "local_orbit_composition.B_Va_1NN-all-combined:" in out
    # etc.


def test_run_fixture_2_collect_event_data(FCCBinaryVacancy_kmc_System, tmp_path):
    """A single run, using a fixture;"""
    system = FCCBinaryVacancy_kmc_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    calculator.collect("selected_event.by_type")
    calculator.collect("selected_event.by_equivalent_index")
    calculator.collect("selected_event.by_equivalent_index_and_direction")
    calculator.collect("selected_event.A_Va_1NN.by_equivalent_index")
    calculator.collect("selected_event.B_Va_1NN.by_equivalent_index")
    calculator.collect("local_orbit_composition.A_Va_1NN-1")
    calculator.collect("local_orbit_composition.B_Va_1NN-1")
    calculator.collect("local_orbit_composition.A_Va_1NN-all-combined")
    calculator.collect("local_orbit_composition.B_Va_1NN-all-combined")
    calculator.collect("local_orbit_composition.A_Va_1NN-all")
    calculator.collect("local_orbit_composition.B_Va_1NN-all")

    # requires re-calculating event states
    calculator.collect("dE_activated.by_type", bin_width=0.01)
    calculator.collect("dE_activated.by_equivalent_index", bin_width=0.01)

    data = calculator.selected_event_function_params.to_dict()
    assert isinstance(data, dict)

    print("SELECTED EVENT FUNCTION PARAMS:")
    print(xtal.pretty_json(data))
    print()

    # TODO: SelectedEventFunctionParams to/from dict parsing test
    selected_event_function_params_in = sampling.SelectedEventFunctionParams.from_dict(
        data
    )
    assert isinstance(
        selected_event_function_params_in, sampling.SelectedEventFunctionParams
    )

    # construct default sampling fixture parameters
    kinetics = calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(output_dir),
        write_observations=True,
    )
    kinetics.sample_by_pass()
    kinetics.sample("clex.formation_energy")
    kinetics.clear_cutoffs()
    kinetics.converge("clex.formation_energy", abs=1e-3)

    print("SAMPLING FIXTURE PARAMS:")
    print(xtal.pretty_json(kinetics.to_dict()))
    print()

    # construct the initial state:
    # start from the default configuration
    # and modify to match mol_composition=[0.899, 0.1, 0.001]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )
    composition_calculator = system.composition_calculator
    assert (
        composition_calculator.num_each_component(
            initial_state.configuration.occupation
        )
        == np.array([899, 100, 1], dtype="int")
    ).all()

    # make_canonical_initial_state should have set the mol_composition
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"], [0.899, 0.1, 0.001]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["param_composition"], [0.1, 0.001]
    )

    # Run
    sampling_fixture = calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=kinetics,
    )
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    # print(sampling_fixture.results.json_samplers["json_step"].values())

    pytest.helpers.validate_summary_file(
        summary_file=summary_file,
        expected_size=1,
        is_canonical=True,
        is_requested_convergence=False,
    )

    # Get the selected event data:
    selected_event_data = calculator.selected_event_data
    assert isinstance(selected_event_data, sampling.SelectedEventData)

    data = selected_event_data.to_dict()
    assert isinstance(data, dict)

    print("SELECTED EVENT DATA:")
    print(xtal.pretty_json(data))
    print()

    expected = [
        "dE_activated.by_equivalent_index",
        "dE_activated.by_type",
        "local_orbit_composition.A_Va_1NN-1",
        "local_orbit_composition.A_Va_1NN-all",
        "local_orbit_composition.A_Va_1NN-all-combined",
        "local_orbit_composition.B_Va_1NN-1",
        "local_orbit_composition.B_Va_1NN-all",
        "local_orbit_composition.B_Va_1NN-all-combined",
        "selected_event.A_Va_1NN.by_equivalent_index",
        "selected_event.B_Va_1NN.by_equivalent_index",
        "selected_event.by_equivalent_index",
        "selected_event.by_equivalent_index_and_direction",
        "selected_event.by_type",
    ]
    histograms_keys = data["histograms"].keys()
    for key in expected:
        assert key in histograms_keys


def test_event_data_summary(FCCBinaryVacancy_kmc_System):
    """Test generating the event data summary"""

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=FCCBinaryVacancy_kmc_System,
    )

    # construct the initial state:
    # start from the default configuration
    # and modify to match mol_composition=[0.899, 0.1, 0.001]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )

    # Set the state and potential
    print("SETTING STATE AND POTENTIAL ... ")
    calculator.set_state_and_potential(state=initial_state)
    print()

    # Set the event data and make the event data summary
    print("SETTING EVENT DATA ... ")
    calculator.make_occ_location()
    calculator.set_event_data()
    event_data_summary = calculator.event_data_summary()
    assert isinstance(event_data_summary, clexmonte.MonteEventDataSummary)
    print()

    # Print the summary
    print("EVENT DATA SUMMARY (PRETTY):")
    print(event_data_summary)
    f = io.StringIO()
    with redirect_stdout(f):
        print(event_data_summary)
    out = f.getvalue()
    assert "Number of unitcells = 1000" in out
    assert "Number of events (total) = 12" in out
    print()

    # Output the summary as a Python dict
    data = event_data_summary.to_dict()
    assert isinstance(data, dict)
    # print("EVENT DATA SUMMARY (JSON):")
    # print(xtal.pretty_json(data))

    expected = [
        "event_list_size",
        "impact_neighborhood",
        "impact_number",
        "mean_time_increment",
        "memory_used",
        "memory_used_MiB",
        "n_events",
        "n_abnormal_events",
        "n_unitcells",
        "rate",
        "stats",
    ]
    keys = list(data.keys())
    for key in expected:
        assert key in keys

    print("DONE\n\n")


def test_SelectedEvent_data(FCCBinaryVacancy_kmc_System_2, tmp_path):
    # A system with formation_energy_eci.2.json chosen to have abnormal events
    system = FCCBinaryVacancy_kmc_System_2
    output_dir = tmp_path / "output"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
        params={},
    )
    calculator.collect("selected_event.by_equivalent_index")
    calculator.collect("dE_activated.by_type", bin_width=0.01)
    calculator.event_data.set_abnormal_event_handling_off()

    # construct default sampling fixture parameters
    kinetics = calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(output_dir),
        write_observations=True,
    )
    kinetics.sample_by_step()
    kinetics.set_max_count(100)
    # kinetics.converge("clex.formation_energy", abs=1e-3)

    selected_events = []

    def print_selected_event_f():
        # Want to print event info? Do something like this:
        # selected_event = calculator.selected_event
        # is_normal = selected_event.event_state.is_normal
        # dE_activated = selected_event.event_state.dE_activated
        # dE_final = selected_event.event_state.dE_final
        # print("~~~")
        # print("dE_activated:", dE_activated)
        # print("dE_final:", dE_final)
        # print("is_normal:", is_normal)
        # print("event_type_name:", event_type_name)
        # print(selected_event.prim_event_data)
        # print(selected_event.event_state)
        # print(selected_event.event_data)
        # print(selected_event)

        selected_event = calculator.selected_event

        # check PrimEventData data:
        data = selected_event.prim_event_data.to_dict()
        assert isinstance(data, dict)

        f = io.StringIO()
        with redirect_stdout(f):
            print(selected_event.prim_event_data)
        out = f.getvalue()
        assert "sites" in out

        # check EventData data:
        data = selected_event.event_data.to_dict()
        assert isinstance(data, dict)

        f = io.StringIO()
        with redirect_stdout(f):
            print(selected_event.event_data)
        out = f.getvalue()
        assert "unitcell_index" in out

        # check EventState data:
        data = selected_event.event_state.to_dict()
        assert isinstance(data, dict)

        f = io.StringIO()
        with redirect_stdout(f):
            print(selected_event.event_state)
        out = f.getvalue()
        assert "dE_activated" in out

        # check SelectedEvent data:
        data = selected_event.to_dict()
        assert isinstance(data, dict)

        f = io.StringIO()
        with redirect_stdout(f):
            print(selected_event)
        out = f.getvalue()
        assert "time_increment" in out

        selected_events.append(data)

    calculator.add_generic_function(
        name="print_event_state",
        description="Print information about the selected event state",
        requires_event_state=False,
        function=print_selected_event_f,
        order=0,
    )
    calculator.evaluate("print_event_state")

    # construct the initial state:
    # start from the default configuration
    # and modify to match mol_composition=[0.899, 0.1, 0.001]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )

    # Run
    calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=kinetics,
    )

    assert len(selected_events) == 100
    for x in selected_events:
        assert isinstance(x, dict)
        assert "event_data" in x
        assert "event_id" in x
        assert "event_index" in x
        assert "event_state" in x
        assert "prim_event_data" in x
        assert "time_increment" in x
        assert "total_rate" in x


def test_write_LocalConfiguration(FCCBinaryVacancy_kmc_System_2, tmp_path):
    """Test writing LocalConfiguration using selected event data"""

    # A system with formation_energy_eci.2.json chosen to have abnormal events
    system = FCCBinaryVacancy_kmc_System_2
    output_dir = tmp_path / "output"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
        params={},
    )
    calculator.collect("selected_event.by_equivalent_index")
    calculator.collect("dE_activated.by_type", bin_width=0.01)
    calculator.event_data.set_abnormal_event_handling_off()

    # construct default sampling fixture parameters
    kinetics = calculator.make_default_sampling_fixture_params(
        label="kinetics",
        output_dir=str(output_dir),
        write_observations=True,
    )
    kinetics.sample_by_step()
    kinetics.set_max_count(100)
    # kinetics.converge("clex.formation_energy", abs=1e-3)

    sampled_local_configs = []

    def print_selected_event_f():
        # Only want barrier-less events? Do something like this:
        # selected_event = calculator.selected_event
        # if selected_event.event_state.is_normal:
        #     return

        # check make_local_configuration:
        x = calculator.make_local_configuration()
        sampled_local_configs.append(x)

        print("# sampled:", len(sampled_local_configs))
        print("local_configuration:")
        print(x)

    calculator.add_generic_function(
        name="print_event_state",
        description="Print information about the selected event state",
        requires_event_state=False,
        function=print_selected_event_f,
        order=0,
    )
    calculator.evaluate("print_event_state")

    # construct the initial state:
    # start from the default configuration
    # and modify to match mol_composition=[0.899, 0.1, 0.001]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 1200.0,
            # one of param/mol composition is needed
            # "param_composition": [0.0, 0.0],
            "mol_composition": [0.899, 0.1, 0.001],
        },
        min_volume=1000,
    )

    # Run
    calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=kinetics,
    )

    # Check the sampled local configurations
    assert len(sampled_local_configs) == 100
    for x in sampled_local_configs:
        assert isinstance(x, LocalConfiguration)
        assert x.configuration is not None
        assert x.pos is not None
        assert x.event_info is not None
