import libcasm.clexmonte as clexmonte
import libcasm.monte as monte


def test_construction_1(Clex_ZrO_Occ_thermo):
    thermo, tmp_path = Clex_ZrO_Occ_thermo
    assert isinstance(thermo, clexmonte.SamplingFixtureParams)

    run_manager = clexmonte.RunManager(
        engine=monte.RandomNumberEngine(),
        sampling_fixture_params=[thermo],
        global_cutoff=True,
    )
    assert isinstance(run_manager, clexmonte.RunManager)
