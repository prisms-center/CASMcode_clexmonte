import json
import os
import os.path
import pathlib
import shutil
import sys

import pytest

import libcasm.xtal as xtal
from libcasm.clexmonte import (
    System,
)
from libcasm.composition import (
    CompositionConverter,
)


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
    original_shared_path = pathlib.Path(os.path.realpath(__file__)).parent / "data"
    session_temp_path = tmpdir_factory.mktemp("session_data")
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
    return CompositionConverter.from_dict(
        FCCBinaryVacancy_system_data["composition_axes"]
    )


@pytest.fixture
def FCCBinaryVacancy_System(FCCBinaryVacancy_system_data, session_shared_datadir):
    return System.from_dict(
        data=FCCBinaryVacancy_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )


@pytest.fixture
def FCCBinaryVacancy_kmc_System(
    FCCBinaryVacancy_kmc_system_data, session_shared_datadir
):
    return System.from_dict(
        data=FCCBinaryVacancy_kmc_system_data,
        search_path=[str(session_shared_datadir / "FCC_binary_vacancy")],
    )
