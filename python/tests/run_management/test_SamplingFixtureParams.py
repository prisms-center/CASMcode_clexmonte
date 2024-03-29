import libcasm.clexmonte as clexmonte


def test_construction_1(Clex_ZrO_Occ_thermo):
    thermo, tmp_path = Clex_ZrO_Occ_thermo
    assert isinstance(thermo, clexmonte.SamplingFixtureParams)
