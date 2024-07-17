import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
import libcasm.xtal as xtal


def test_constructors_1(Clex_ZrO_Occ_System):
    system = Clex_ZrO_Occ_System

    # The mol composition element meaning is determined by the
    # order of components in the composition calculator
    assert system.composition_calculator.components() == ["Zr", "Va", "O"]

    # The param_composition meaning is determined by the origin and end member
    # mol compositions
    assert np.allclose(system.composition_converter.origin(), [2.0, 2.0, 0.0])
    assert np.allclose(system.composition_converter.end_member(0), [2.0, 0.0, 2.0])

    calculator = clexmonte.MonteCalculator(
        method="canonical",
        system=system,
    )
    assert isinstance(calculator, clexmonte.MonteCalculator)
    with pytest.raises(Exception):
        assert isinstance(calculator.potential, clexmonte.MontePotential)
    with pytest.raises(Exception):
        assert isinstance(calculator.state_data, clexmonte.StateData)

    # default configuration is occupied by Va: [2.0, 2.0, 0.0], which corresponds
    # to the origin composition as defined in the system's composition axes
    state = clexmonte.MonteCarloState(
        configuration=system.make_default_configuration(
            transformation_matrix_to_super=np.eye(3, dtype="int") * 2,
        ),
        conditions={
            "temperature": 300.0,
            "param_composition": [0.0],  # <- one of param/mol composition is needed
            # "mol_composition": [2.0, 2.0, 0.0],
        },
    )
    composition_calculator = system.composition_calculator
    composition_converter = system.composition_converter
    mol_composition = composition_calculator.mean_num_each_component(
        state.configuration.occupation
    )
    assert np.allclose(mol_composition, [2.0, 2.0, 0.0])
    # assert np.allclose(
    #     mol_composition, state.conditions.vector_values["mol_composition"]
    # )
    param_composition = composition_converter.param_composition(mol_composition)
    assert np.allclose(param_composition, [0.0])
    assert np.allclose(
        param_composition, state.conditions.vector_values["param_composition"]
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
        method="canonical",
        system=system,
    )

    # construct default sampling fixture parameters
    thermo = calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )
    print(xtal.pretty_json(thermo.to_dict()))

    # construct the initial state:
    # start from the default configuration and modify to match param_composition=[0.5]
    initial_state, motif = clexmonte.make_canonical_initial_state(
        calculator=calculator,
        conditions={
            "temperature": 300.0,
            "param_composition": [0.5],  # <- one of param/mol composition is needed
            # "mol_composition": [2.0, 1.0, 1.0],
        },
        min_volume=1000,
    )
    composition_calculator = system.composition_calculator
    assert np.allclose(
        composition_calculator.mean_num_each_component(
            initial_state.configuration.occupation
        ),
        [2.0, 1.0, 1.0],
    )
    # make_canonical_initial_state should have set the mol_composition
    assert np.allclose(
        initial_state.conditions.vector_values["mol_composition"], [2.0, 1.0, 1.0]
    )
    assert np.allclose(
        initial_state.conditions.vector_values["param_composition"], [0.5]
    )

    # Run
    sampling_fixture = calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=thermo,
    )
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    pytest.helpers.validate_summary_file(
        summary_file=summary_file,
        expected_size=1,
        is_canonical=True,
    )
