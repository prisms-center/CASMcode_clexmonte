import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
import libcasm.monte.events as monte_events
import libcasm.monte.sampling as sampling
import libcasm.xtal as xtal


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
    state = clexmonte.MonteCarloState(
        configuration=system.make_default_configuration(
            transformation_matrix_to_super=np.eye(3, dtype="int") * 2,
        ),
        conditions={
            "temperature": 300.0,
            "param_composition": [0.0, 0.0],  # <-one of param/mol composition is needed
            # "mol_composition": [1.0, 0.0, 0.0],
        },
    )
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
    calculator.set_state_and_potential(state=state)
    assert isinstance(calculator.potential, clexmonte.MontePotential)
    assert isinstance(calculator.state_data, clexmonte.StateData)
    assert isinstance(calculator.event_data, clexmonte.MonteEventData)

    ## StateData constructor
    state_data = clexmonte.StateData(
        system=system,
        state=state,
        occ_location=None,
    )
    assert isinstance(state_data, clexmonte.StateData)

    ## MontePotential constructor
    potential = clexmonte.MontePotential(calculator=calculator, state=state)
    assert isinstance(potential, clexmonte.MontePotential)

    # ## Set event data
    with pytest.raises(Exception):
        # no occ_location set
        calculator.set_event_data()

    occ_location = calculator.make_occ_location()
    assert isinstance(occ_location, monte_events.OccLocation)
    assert isinstance(calculator.event_data, clexmonte.MonteEventData)

    ## MonteEventData constructor
    event_data = clexmonte.MonteEventData(calculator=calculator, state=state)
    assert isinstance(event_data, clexmonte.MonteEventData)


def test_run_fixture_1(FCCBinaryVacancy_kmc_System, tmp_path):
    """A single run, using a fixture"""
    system = FCCBinaryVacancy_kmc_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a semi-grand canonical MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method="kinetic",
        system=system,
    )

    def print_step_f():
        step = calculator.kinetics_data.sampling_fixture.step
        # print("step:", step)
        return np.array([step])

    print_step_sampling_f = sampling.StateSamplingFunction(
        name="print_step",
        description="test",
        shape=[],
        function=print_step_f,
        component_names=["step"],
    )

    def json_step_f():
        json_step = {"step": calculator.kinetics_data.sampling_fixture.step}
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
