"""Test C++ implemented property ising_cpp with Python implemented Monte Carlo loop"""
import json
import numpy as np

import libcasm.clexmonte as clexmonte
import libcasm.clexmonte.semigrand_canonical as sgc
import libcasm.monte as monte
import libcasm.monte.sampling as monte_sampling
import libcasm.xtal as xtal


def validate_summary_data(subdata: dict, expected_keys: list[str], expected_size: int):
    for x in expected_keys:
        assert x in subdata
        if "component_names" in subdata[x]:
            # non-scalar analysis functions & conditions
            for y in subdata[x]["component_names"]:
                assert len(subdata[x][y]) == expected_size
        elif "value" in subdata[x]:
            # scalar analysis functions
            assert subdata[x]["shape"] == []
            assert len(subdata[x]["value"]) == expected_size
        else:
            # completion_check_params
            assert len(subdata[x]) == expected_size


def validate_statistics_data(
    subdata: dict,
    expected_keys: list[str],
    expected_size: int,
):
    for x in expected_keys:
        assert x in subdata
        if "component_names" in subdata[x]:
            for y in subdata[x]["component_names"]:
                assert y in subdata[x]
                for z in ["mean", "calculated_precision"]:
                    assert z in subdata[x][y]
                    assert len(subdata[x][y][z]) == expected_size
        else:
            assert subdata[x]["shape"] == []
            assert "value" in subdata[x]
            for z in ["mean", "calculated_precision"]:
                assert z in subdata[x]["value"]
                assert len(subdata[x]["value"][z]) == expected_size


def test_constructors_1(Clex_ZrO_Occ_System):
    system = Clex_ZrO_Occ_System

    potential = sgc.SemiGrandCanonicalPotential(system=system)
    assert isinstance(potential, sgc.SemiGrandCanonicalPotential)

    conditions = sgc.SemiGrandCanonicalConditions(
        composition_converter=system.composition_converter
    )
    assert isinstance(conditions, sgc.SemiGrandCanonicalConditions)

    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)
    assert isinstance(mc_calculator, sgc.SemiGrandCanonicalCalculator)


def test_run_fixture_0(Clex_ZrO_Occ_System, tmp_path):
    system = Clex_ZrO_Occ_System
    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a SemiGrandCanonicalCalculator
    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)

    # construct default sampling fixture parameters
    thermo = mc_calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )
    print(xtal.pretty_json(thermo.to_dict()))

    # construct the initial state (default configuration)
    initial_state, motif, motif_id = mc_calculator.make_initial_state(
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [-1.0],
        },
        min_volume=1000,
    )

    # Run
    sampling_fixture = mc_calculator.run_fixture(
        state=initial_state,
        sampling_fixture_params=thermo,
    )
    assert isinstance(sampling_fixture, clexmonte.SamplingFixture)
    assert summary_file.exists() and summary_file.is_file()
    with open(summary_file, "r") as f:
        data = json.load(f)
    print(xtal.pretty_json(data))

    expected_size = 1
    assert "analysis" in data
    validate_summary_data(
        subdata=data["analysis"],
        expected_keys=[
            "heat_capacity",
            "mol_susc",
            "param_susc",
            "mol_thermochem_susc",
            "param_thermochem_susc",
        ],
        expected_size=expected_size,
    )

    assert "completion_check_results" in data
    validate_summary_data(
        subdata=data["completion_check_results"],
        expected_keys=[
            "N_samples",
            "N_samples_for_all_to_equilibrate",
            "N_samples_for_statistics",
            "acceptance_rate",
            "all_converged",
            "all_equilibrated",
            "count",
            "elapsed_clocktime",
        ],
        expected_size=expected_size,
    )

    assert "conditions" in data
    validate_summary_data(
        subdata=data["conditions"],
        expected_keys=["temperature", "param_chem_pot"],
        expected_size=expected_size,
    )

    assert "statistics" in data
    validate_statistics_data(
        subdata=data["statistics"],
        expected_keys=[
            "formation_energy",
            "mol_composition",
            "param_composition",
            "potential_energy",
        ],
        expected_size=expected_size,
    )
    assert "is_converged" in data["statistics"]["potential_energy"]["value"]
    assert "is_converged" in data["statistics"]["param_composition"]["a"]


def test_run_fixture_1(Clex_ZrO_Occ_System, tmp_path):
    system = Clex_ZrO_Occ_System

    # path = session_shared_datadir / "Clex_ZrO_Occ" / "system.json"
    # with open(path, "r") as f:
    #     Clex_ZrO_Occ_system_data = json.load(f)
    #
    # system = System.from_dict(
    #     data=Clex_ZrO_Occ_system_data,
    #     search_path=[str(session_shared_datadir / "Clex_ZrO_Occ")],
    # )

    output_dir = tmp_path / "output"
    summary_file = output_dir / "summary.json"

    # construct a SemiGrandCanonicalCalculator
    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)

    # construct default sampling fixture parameters
    thermo = mc_calculator.make_default_sampling_fixture_params(
        label="thermo",
        output_dir=str(output_dir),
    )

    # set lower convergence level for potential_energy
    thermo.converge(quantity="potential_energy", abs=2e-3)

    # set lower convergence level for param_composition("a")
    thermo.converge(quantity="param_composition", abs=2e-3, component_name=["a"])

    # construct the initial state (default configuration)
    state, motif, motif_id = mc_calculator.make_initial_state(
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
        sampling_fixture = mc_calculator.run_fixture(
            state=state,
            sampling_fixture_params=thermo,
        )
        assert isinstance(sampling_fixture, clexmonte.SamplingFixture)

    assert summary_file.exists() and summary_file.is_file()
    with open(summary_file, "r") as f:
        data = json.load(f)
    print(xtal.pretty_json(data))

    expected_size = len(x_list)
    assert "analysis" in data
    validate_summary_data(
        subdata=data["analysis"],
        expected_keys=[
            "heat_capacity",
            "mol_susc",
            "param_susc",
            "mol_thermochem_susc",
            "param_thermochem_susc",
        ],
        expected_size=expected_size,
    )

    assert "completion_check_results" in data
    validate_summary_data(
        subdata=data["completion_check_results"],
        expected_keys=[
            "N_samples",
            "N_samples_for_all_to_equilibrate",
            "N_samples_for_statistics",
            "acceptance_rate",
            "all_converged",
            "all_equilibrated",
            "count",
            "elapsed_clocktime",
        ],
        expected_size=expected_size,
    )

    assert "conditions" in data
    validate_summary_data(
        subdata=data["conditions"],
        expected_keys=["temperature", "param_chem_pot"],
        expected_size=expected_size,
    )

    assert "statistics" in data
    validate_statistics_data(
        subdata=data["statistics"],
        expected_keys=[
            "formation_energy",
            "mol_composition",
            "param_composition",
            "potential_energy",
        ],
        expected_size=expected_size,
    )
    assert "is_converged" in data["statistics"]["potential_energy"]["value"]
    assert "is_converged" in data["statistics"]["param_composition"]["a"]


def test_run_1(Clex_ZrO_Occ_System, tmp_path):
    system = Clex_ZrO_Occ_System

    # construct a SemiGrandCanonicalCalculator
    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)

    assert isinstance(mc_calculator, sgc.SemiGrandCanonicalCalculator)

    # construct sampling functions
    sampling_functions = mc_calculator.standard_sampling_functions()
    json_sampling_functions = mc_calculator.standard_json_sampling_functions()
    analysis_functions = mc_calculator.standard_analysis_functions()

    assert isinstance(sampling_functions, monte_sampling.StateSamplingFunctionMap)
    assert isinstance(
        json_sampling_functions, monte_sampling.jsonStateSamplingFunctionMap
    )

    # construct the initial state
    initial_state = system.make_default_state(
        transformation_matrix_to_super=np.array(
            [
                [10, 0, 0],
                [0, 10, 0],
                [0, 0, 10],
            ]
        ),
    )
    initial_state.conditions.scalar_values["temperature"] = 300.0
    initial_state.conditions.vector_values["param_chem_pot"] = [0.0]

    assert isinstance(initial_state, clexmonte.MonteCarloState)

    ### thermo sampling fixture ###

    # completion check params
    completion_check_params = monte_sampling.CompletionCheckParams()
    completion_check_params.cutoff_params.min_sample = 100
    completion_check_params.log_spacing = False
    completion_check_params.check_begin = 100
    completion_check_params.check_period = 10

    # Set requested precision
    monte_sampling.converge(
        sampling_functions,
        completion_check_params,
    ).set_precision(
        "potential_energy",
        abs=0.001,
    ).set_precision(
        "param_composition",
        abs=0.001,
    )

    sampling_params = monte_sampling.SamplingParams(
        sampler_names=[
            "mol_composition",
            "param_composition",
            "potential_energy",
        ],
    )
    thermo_sampling_fixture_params = clexmonte.SamplingFixtureParams(
        label="thermo",
        sampling_functions=sampling_functions,
        json_sampling_functions=json_sampling_functions,
        analysis_functions=analysis_functions,
        sampling_params=sampling_params,
        completion_check_params=completion_check_params,
        output_dir=str(tmp_path / "output"),
    )
    ###

    run_manager = clexmonte.RunManager(
        engine=monte.RandomNumberEngine(),
        sampling_fixture_params=[thermo_sampling_fixture_params],
        global_cutoff=True,
    )

    # Run
    mc_calculator.run(
        state=initial_state,
        run_manager=run_manager,
    )

    assert isinstance(run_manager, clexmonte.RunManager)


def test_run_2(Clex_ZrO_Occ_System, tmp_path):
    system = Clex_ZrO_Occ_System

    # construct a SemiGrandCanonicalCalculator
    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)

    # construct default sampling fixture parameters
    thermo = mc_calculator.make_default_sampling_fixture_params(
        label="thermo",
        write_observations=True,
        write_trajectory=True,
        output_dir=str(tmp_path / "output"),
    )

    run_manager = clexmonte.RunManager(
        engine=monte.RandomNumberEngine(),
        sampling_fixture_params=[thermo],
        global_cutoff=True,
    )

    # construct the initial state (default configuration)
    initial_state = system.make_default_state(
        transformation_matrix_to_super=np.array(
            [
                [10, 0, 0],
                [0, 10, 0],
                [0, 0, 10],
            ]
        ),
    )
    initial_state.conditions.scalar_values["temperature"] = 300.0
    initial_state.conditions.vector_values["param_chem_pot"] = [0.0]

    assert isinstance(initial_state, clexmonte.MonteCarloState)

    # Run
    mc_calculator.run(
        state=initial_state,
        run_manager=run_manager,
    )

    assert isinstance(run_manager, clexmonte.RunManager)


def test_run_3(Clex_ZrO_Occ_System, tmp_path):
    system = Clex_ZrO_Occ_System

    # construct a SemiGrandCanonicalCalculator
    mc_calculator = sgc.SemiGrandCanonicalCalculator(system=system)

    # construct default sampling fixture parameters
    thermo = mc_calculator.make_default_sampling_fixture_params(
        label="thermo",
        write_observations=True,
        write_trajectory=True,
        output_dir=str(tmp_path / "output"),
    )

    run_manager = clexmonte.RunManager(
        engine=monte.RandomNumberEngine(),
        sampling_fixture_params=[thermo],
        global_cutoff=True,
    )

    # construct the initial state (default configuration)
    initial_state, motif, id = sgc.make_initial_state(
        system=system,
        conditions={
            "temperature": 300.0,
            "param_chem_pot": [0.0],
        },
        # dirs = "abc",
        min_volume=1000,
        transformation_matrix_to_super=None,
        motif=None,
        configurations=None,
    )
    assert isinstance(initial_state, clexmonte.MonteCarloState)

    # Run
    mc_calculator.run(
        state=initial_state,
        run_manager=run_manager,
    )

    assert isinstance(run_manager, clexmonte.RunManager)
