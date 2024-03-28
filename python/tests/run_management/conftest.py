import numpy as np
import pytest

import libcasm.clexmonte as clexmonte
import libcasm.clexmonte.semigrand_canonical as sgc
import libcasm.monte.sampling as monte_sampling


@pytest.fixture
def Clex_ZrO_Occ_thermo(Clex_ZrO_Occ_System, tmp_path):
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
    thermo = clexmonte.SamplingFixtureParams(
        label="thermo",
        sampling_functions=sampling_functions,
        json_sampling_functions=json_sampling_functions,
        analysis_functions=analysis_functions,
        sampling_params=sampling_params,
        completion_check_params=completion_check_params,
        output_dir=str(tmp_path / "output"),
    )

    return (thermo, tmp_path)