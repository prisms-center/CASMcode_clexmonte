import pathlib
import typing

import libcasm.clexmonte as clexmonte
import libcasm.monte.sampling as sampling
from libcasm.clexmonte.site_iterators import (
    AsymmetricUnitSiteIterator,
    ChemicalTypeSiteIterator,
    CompleteSiteIterator,
    OccCandidateSiteIterator,
    OccupantTypeSiteIterator,
    SublatticeSiteIterator,
)


def add_sampler(
    calculator: clexmonte.MonteCalculator,
    sampling_fixture_params: clexmonte.SamplingFixture,
    cls: typing.Any,
    kwargs: dict,
):
    site_iterator = cls(
        calculator=calculator,
        **kwargs,
    )

    def _function():
        # configuration = calculator.state_data.state.configuration

        # occ_location = calculator.state_data.occ_location
        # print()
        # for cand in occ_location.candidate_list():
        #     a = cand.asymmetric_unit_index
        #     s = cand.species_index
        #     print(f"cand: (a={a}, s={s}")
        # print()
        # exit()

        n_visited = 0
        for site in site_iterator:
            # Any Monte Carlo:
            # print("site.l:", site.l, "occ: ", configuration.occ(site.l))
            # print("  site.mol_id:", site.mol_id)

            # KMC w/ atomic occupants:
            # print("  site.atom_id:", site.atom_id)
            # print("  site.atom_n_jumps:", site.atom_n_jumps)
            # print("  site.atom_coordinate_frac:", site.atom_coordinate_frac)
            # print("  site.atom_coordinate_cart:", site.atom_coordinate_cart)
            n_visited += 1

        # print("done")
        # print()
        # sys.stdout.flush()

        return [n_visited]

    sampling_f = sampling.StateSamplingFunction(
        name="n_visited",
        description="desc",
        shape=[],
        function=_function,
        component_names=["n_visited"],
    )

    sampling_fixture_params.add_sampling_function(sampling_f)
    sampling_fixture_params.sample("n_visited")


def run(
    system: clexmonte.System,
    cls: typing.Any,
    kwargs: dict,
    tmp_path: pathlib.Path,
    method: str,
    conditions: dict,
):
    """Run a MonteCalculator and use a site iterator in a custom sampling function

    Parameters
    ----------
    system: clexmonte.System
        The system.
    cls: Any
        The site iterator class
    kwargs:
        The site iterator class constructor arguments, excluding `calculator`.
    tmp_path: pathlib.Path
        A temporary directory for output.
    method: str
        The Monte Carlo method to use ("canonical" or "semigrand_canonical").

    Returns
    -------
    sampling_fixture: clexmonte.SamplingFixture
        The sampling fixture with results.
    """

    output_dir = tmp_path / "output"

    # construct a MonteCalculator
    calculator = clexmonte.MonteCalculator(
        method=method,
        system=system,
        params={"verbosity": "quiet"} if method == "kinetic" else {},
    )

    # construct default sampling fixture parameters
    fixture_params = calculator.make_default_sampling_fixture_params(
        label="data",
        output_dir=str(output_dir),
    )
    fixture_params.set_min_count(100)
    fixture_params.set_max_count(100)
    fixture_params.sample_by_pass(begin=0.0, period=1.0)
    # print(pretty_json(fixture_params.to_dict()))

    # construct the initial state (default configuration)
    vol = 1000
    if method in ["canonical", "kinetic"]:
        state, motif = clexmonte.make_canonical_initial_state(
            calculator=calculator,
            conditions=conditions,
            min_volume=vol,
        )
    elif method == "semigrand_canonical":
        state, motif, motif_id = clexmonte.make_initial_state(
            calculator=calculator,
            conditions=conditions,
            min_volume=vol,
        )
    else:
        raise ValueError(f"Unknown method: {method}")
    # print(pretty_json(state.to_dict()))

    add_sampler(
        calculator=calculator,
        sampling_fixture_params=fixture_params,
        cls=cls,
        kwargs=kwargs,
    )

    # Run and return sampling fixture
    return calculator.run_fixture(
        state=state,
        sampling_fixture_params=fixture_params,
    )


def test_canonical_site_iterators_1(FCCBinaryVacancy_System, tmpdir):
    system = FCCBinaryVacancy_System

    Va_comp = 4.0 / 1000.0
    B_comp = 0.2
    A_comp = 1.0 - B_comp - Va_comp
    conditions = {
        "temperature": 300.0,
        "mol_composition": [A_comp, B_comp, Va_comp],
    }

    # 1: CompleteSiteIterator
    sampling_fixture = run(
        system=system,
        cls=CompleteSiteIterator,
        kwargs={},
        tmp_path=tmpdir.mkdir("complete"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 2: SublatticeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=SublatticeSiteIterator,
        kwargs={"sublattice": [0]},
        tmp_path=tmpdir.mkdir("sublattice"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 3: AsymmetricUnitSiteIterator
    sampling_fixture = run(
        system=system,
        cls=AsymmetricUnitSiteIterator,
        kwargs={"asym": [0]},
        tmp_path=tmpdir.mkdir("asymmetric_unit"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 4: OccCandidateSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccCandidateSiteIterator,
        kwargs={"candidate_indices": [2]},
        tmp_path=tmpdir.mkdir("occ_candidate"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)

    # 5: ChemicalTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=ChemicalTypeSiteIterator,
        kwargs={"chemical_name": ["Va"]},
        tmp_path=tmpdir.mkdir("chemical_type"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)

    # 6: OccupantTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccupantTypeSiteIterator,
        kwargs={"occupant_name": ["Va"]},
        tmp_path=tmpdir.mkdir("occupant_type"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)


def test_canonical_site_iterators_2(Clex_ZrO_Occ_System, tmpdir):
    system = Clex_ZrO_Occ_System

    Zr_comp = 2.0
    Va_comp = 1667.0 / 1000.0
    O_comp = 2.0 - Va_comp

    conditions = {
        "temperature": 300.0,
        "mol_composition": [Zr_comp, Va_comp, O_comp],
    }

    # 1: CompleteSiteIterator
    sampling_fixture = run(
        system=system,
        cls=CompleteSiteIterator,
        kwargs={},
        tmp_path=tmpdir.mkdir("complete"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4 * 1000.0)

    # 2: SublatticeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=SublatticeSiteIterator,
        kwargs={"sublattice": [2, 3]},
        tmp_path=tmpdir.mkdir("sublattice"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 2 * 1000.0)

    # 3a: AsymmetricUnitSiteIterator
    sampling_fixture = run(
        system=system,
        cls=AsymmetricUnitSiteIterator,
        kwargs={"asym": [1]},
        tmp_path=tmpdir.mkdir("asymmetric_unit_a"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 2 * 1000.0)

    # 3b: AsymmetricUnitSiteIterator
    sampling_fixture = run(
        system=system,
        cls=AsymmetricUnitSiteIterator,
        kwargs={"asym": [0, 1]},
        tmp_path=tmpdir.mkdir("asymmetric_unit_b"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4 * 1000.0)

    # 4a: OccCandidateSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccCandidateSiteIterator,
        kwargs={"candidate_indices": [0]},
        tmp_path=tmpdir.mkdir("occ_candidate_a"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1667.0)

    # 4b: OccCandidateSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccCandidateSiteIterator,
        kwargs={"candidate_indices": [0, 1]},
        tmp_path=tmpdir.mkdir("occ_candidate_b"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 2000.0)

    # 5: ChemicalTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=ChemicalTypeSiteIterator,
        kwargs={"chemical_name": ["Va"]},
        tmp_path=tmpdir.mkdir("chemical_type"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1667.0)

    # 6: OccupantTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccupantTypeSiteIterator,
        kwargs={"occupant_name": ["Va"]},
        tmp_path=tmpdir.mkdir("occupant_type"),
        method="canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1667.0)


def test_semigrand_canonical_site_iterators(FCCBinaryVacancy_System, tmpdir):
    system = FCCBinaryVacancy_System

    conditions = {
        "temperature": 300.0,
        "param_chem_pot": [0.0, 0.22],
    }

    # 1: CompleteSiteIterator
    sampling_fixture = run(
        system=system,
        cls=CompleteSiteIterator,
        kwargs={},
        tmp_path=tmpdir.mkdir("complete"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 2: SublatticeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=SublatticeSiteIterator,
        kwargs={"sublattice": [0]},
        tmp_path=tmpdir.mkdir("sublattice"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 3: AsymmetricUnitSiteIterator
    sampling_fixture = run(
        system=system,
        cls=AsymmetricUnitSiteIterator,
        kwargs={"asym": [0]},
        tmp_path=tmpdir.mkdir("asymmetric_unit"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 4: OccCandidateSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccCandidateSiteIterator,
        kwargs={"candidate_indices": [2]},
        tmp_path=tmpdir.mkdir("occ_candidate"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x >= 0.0) and all(x <= 1000.0)

    # 5: ChemicalTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=ChemicalTypeSiteIterator,
        kwargs={"chemical_name": ["Va"]},
        tmp_path=tmpdir.mkdir("chemical_type"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x >= 0.0) and all(x <= 1000.0)

    # 6: OccupantTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccupantTypeSiteIterator,
        kwargs={"occupant_name": ["Va"]},
        tmp_path=tmpdir.mkdir("occupant_type"),
        method="semigrand_canonical",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x >= 0.0) and all(x <= 1000.0)


def test_kinetic_site_iterators(FCCBinaryVacancy_kmc_System, tmpdir):
    system = FCCBinaryVacancy_kmc_System

    Va_comp = 4.0 / 1000.0
    B_comp = 0.2
    A_comp = 1.0 - B_comp - Va_comp
    conditions = {
        "temperature": 300.0,
        "mol_composition": [A_comp, B_comp, Va_comp],
    }

    # 1: CompleteSiteIterator
    sampling_fixture = run(
        system=system,
        cls=CompleteSiteIterator,
        kwargs={},
        tmp_path=tmpdir.mkdir("complete"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 2: SublatticeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=SublatticeSiteIterator,
        kwargs={"sublattice": [0]},
        tmp_path=tmpdir.mkdir("sublattice"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 3: AsymmetricUnitSiteIterator
    sampling_fixture = run(
        system=system,
        cls=AsymmetricUnitSiteIterator,
        kwargs={"asym": [0]},
        tmp_path=tmpdir.mkdir("asymmetric_unit"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 1000.0)

    # 4: OccCandidateSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccCandidateSiteIterator,
        kwargs={"candidate_indices": [2]},
        tmp_path=tmpdir.mkdir("occ_candidate"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)

    # 5: ChemicalTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=ChemicalTypeSiteIterator,
        kwargs={"chemical_name": ["Va"]},
        tmp_path=tmpdir.mkdir("chemical_type"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)

    # 6: OccupantTypeSiteIterator
    sampling_fixture = run(
        system=system,
        cls=OccupantTypeSiteIterator,
        kwargs={"occupant_name": ["Va"]},
        tmp_path=tmpdir.mkdir("occupant_type"),
        method="kinetic",
        conditions=conditions,
    )

    x = sampling_fixture.results.samplers["n_visited"].values()

    # assert that x is a shape (101,1) array with all elements == 4.0
    assert x.shape == (101, 1)
    assert all(x == 4.0)
