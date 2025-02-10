import json
import os
import os.path
import pathlib
import shutil
import sys

import pytest

import libcasm.clexmonte as clexmonte
import libcasm.composition as casmcomp
import libcasm.xtal as xtal


def _win32_longpath(path):
    """
    Helper function to add the long path prefix for Windows, so that shutil.copytree
     won't fail while working with paths with 255+ chars.
    """
    if sys.platform == "win32":
        # The use of os.path.normpath here is necessary since "the "\\?\" prefix
        # to a path string tells the Windows APIs to disable all string parsing
        # and to send the string that follows it straight to the file system".
        # (See https://docs.microsoft.com/pt-br/windows/desktop/FileIO/naming-a-file)
        return "\\\\?\\" + os.path.normpath(path)
    else:
        return path


@pytest.fixture(scope="session")
def session_shared_datadir(tmpdir_factory):
    tests_path = pathlib.Path(os.path.realpath(__file__)).parent
    original_shared_path = tests_path / "data"
    # session_temp_path = tmpdir_factory.mktemp("session_data")
    session_temp_path = tests_path / "data.tmp"
    shutil.copytree(
        _win32_longpath(original_shared_path),
        _win32_longpath(str(session_temp_path)),
        dirs_exist_ok=True,
    )
    return session_temp_path


@pytest.fixture
def FCCBinaryVacancy_system_data(session_shared_datadir):
    path = session_shared_datadir / "FCC_binary_vacancy" / "system.json"
    with open(path, "r") as f:
        return json.load(f)


@pytest.fixture
def FCCBinaryVacancy_kmc_system_data(session_shared_datadir):
    path = session_shared_datadir / "FCC_binary_vacancy" / "kmc_system.json"
    with open(path, "r") as f:
        return json.load(f)


@pytest.fixture
def FCCBinaryVacancy_xtal_prim(FCCBinaryVacancy_system_data):
    return xtal.Prim.from_dict(FCCBinaryVacancy_system_data["prim"])


@pytest.fixture
def FCCBinaryVacancy_CompositionConverter(FCCBinaryVacancy_system_data):
    return casmcomp.CompositionConverter.from_dict(
        FCCBinaryVacancy_system_data["composition_axes"]
    )


@pytest.fixture
def FCCBinaryVacancy_System(FCCBinaryVacancy_system_data, session_shared_datadir):
    return clexmonte.System.from_dict(
        data=FCCBinaryVacancy_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )


@pytest.fixture
def FCCBinaryVacancy_kmc_System(
    FCCBinaryVacancy_kmc_system_data, session_shared_datadir
):
    return clexmonte.System.from_dict(
        data=FCCBinaryVacancy_kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )


@pytest.fixture
def FCCBinaryVacancy_kmc_System_2(session_shared_datadir):
    """A KMC system with formation_energy_eci.2.json chosen to have events with no
    barriers
    """
    path = session_shared_datadir / "FCC_binary_vacancy" / "kmc_system.2.json"
    with open(path, "r") as f:
        kmc_system_data = json.load(f)
    return clexmonte.System.from_dict(
        data=kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )


@pytest.fixture
def Clex_ZrO_Occ_system_data(session_shared_datadir):
    path = session_shared_datadir / "Clex_ZrO_Occ" / "system.json"
    with open(path, "r") as f:
        return json.load(f)


# function to add two numbers


@pytest.fixture
def Clex_ZrO_Occ_System(Clex_ZrO_Occ_system_data, session_shared_datadir):
    return clexmonte.System.from_dict(
        data=Clex_ZrO_Occ_system_data,
        search_path=[str(session_shared_datadir / "Clex_ZrO_Occ")],
    )


@pytest.helpers.register
class CalculatorTestRunner:
    """A class to help with testing MonteCalculators."""

    def __init__(
        self,
        system: clexmonte.System,
        method: str,
        params: dict,
        output_dir: pathlib.Path,
    ):
        self.system = system
        self.output_dir = output_dir
        self.calculator = clexmonte.MonteCalculator(
            method=method,
            system=system,
            params=params,
        )


@pytest.helpers.register
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


@pytest.helpers.register
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


@pytest.helpers.register
def validate_summary_file(
    summary_file: pathlib.Path,
    expected_size: int,
    is_canonical: bool = False,
    is_requested_convergence: bool = True,
):
    assert summary_file.exists() and summary_file.is_file()
    with open(summary_file, "r") as f:
        data = json.load(f)
    print(xtal.pretty_json(data))

    if is_canonical is False:
        expected_conditions_keys = ["temperature", "param_chem_pot"]
        expected_analysis_keys = [
            "heat_capacity",
            "mol_susc",
            "param_susc",
            "mol_thermochem_susc",
            "param_thermochem_susc",
        ]
    else:
        expected_conditions_keys = [
            "temperature",
            "param_composition",
            "mol_composition",
        ]
        expected_analysis_keys = ["heat_capacity"]

    assert "analysis" in data
    validate_summary_data(
        subdata=data["analysis"],
        expected_keys=expected_analysis_keys,
        expected_size=expected_size,
    )

    assert "completion_check_results" in data
    expected_completion_check_results_keys = [
        "N_samples",
        "N_samples_for_statistics",
        "acceptance_rate",
        "count",
        "elapsed_clocktime",
    ]
    if is_requested_convergence is True:
        expected_completion_check_results_keys += [
            "N_samples_for_all_to_equilibrate",
            "all_converged",
            "all_equilibrated",
        ]
    validate_summary_data(
        subdata=data["completion_check_results"],
        expected_keys=expected_completion_check_results_keys,
        expected_size=expected_size,
    )

    assert "conditions" in data
    validate_summary_data(
        subdata=data["conditions"],
        expected_keys=expected_conditions_keys,
        expected_size=expected_size,
    )

    assert "statistics" in data
    validate_statistics_data(
        subdata=data["statistics"],
        expected_keys=[
            "clex.formation_energy",
            "mol_composition",
            "param_composition",
            "potential_energy",
        ],
        expected_size=expected_size,
    )

    if is_requested_convergence is True:
        if is_canonical is False:
            assert "is_converged" in data["statistics"]["potential_energy"]["value"]
            assert "is_converged" in data["statistics"]["param_composition"]["a"]
        else:
            assert "is_converged" in data["statistics"]["potential_energy"]["value"]
