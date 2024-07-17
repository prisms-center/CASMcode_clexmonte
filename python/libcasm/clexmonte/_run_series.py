# """Run a series of Monte Carlo simulations"""
# import json
# from typing import Union
# import pathlib
#
# import libcasm.monte as monte
# import libcasm.xtal as xtal
#
#
# from ._clexmonte_monte_calculator import MonteCalculator
# from ._clexmonte_run_management import RunManager
#
# from ._FixedConfigGenerator import FixedConfigGenerator
# from ._IncrementalConditionsStateGenerator import IncrementalConditionsStateGenerator
# from ._RunData import (
#     RunData,
#     RunDataOutputParams,
# )
# from libcasm.monte import MethodLog
# from libcasm.monte.events import OccLocation
# from libcasm.configuration import Configuration, Supercell
#
#
# def run_series_a(
#     system_path: Union[str, pathlib.Path],
#     run_params_path: Union[str, pathlib.Path],
#     search_path: list[Union[str, pathlib.Path]],
# ):
#     """Run a series of Monte Carlo simulations
#
#     Parameters
#     ----------
#     system_path: Union[str, pathlib.Path]
#         The path to a JSON formatted file specifying the
#         :class:`~libcasm.clexmonte.System`.
#     run_params_path: Union[str, pathlib.Path]
#         The path to a JSON formatted file specifying the Monte Carlo simulation run
#         parameters.
#     :param search_path:
#     :return:
#     """
#     with open(system_path, "r") as f:
#         system_data = json.load(f)
#     system = System.from_dict(
#         data=system_data,
#         search_path=search_path,
#     )
#
#
# def run_series(
#     calculator: MonteCalculator,
#     supercell: Supercell,
#     motif: Configuration,
#     run_manager: RunManager,
#     initial_conditions: Union[dict, monte.ValueMap],
#     conditions_increment: Union[dict, monte.ValueMap],
#     n_states: int,
#     output_dir: Union[str, pathlib.Path],
#     dependent_runs: bool,
# ) -> RunManager:
#     """Perform a series of Monte Carlo calculations, evolving the input state
#
#     Parameters
#     ----------
#     calculator : MonteCalculator
#         The Monte Carlo calculator.
#     supercell: libcasm.configuration.Supercell
#         The Monte Carlo simulation supercell.
#     motif: libcasm.configuration.Configuration
#         Initial configuration, which will be copied and tiled into the
#         Monte Carlo supercell. If a perfect tiling can be made by
#         applying factor group operations, a note is printed indicating
#         which operation is applied. A warning is printed if there is
#         no perfect tiling and the `motif` is used without reorientation
#         to fill the supercell imperfectly. If `supercell` is given but
#         no `motif` is provided, the default configuration is used.
#     run_manager: RunManager
#         Specifies sampling and convergence criteria and collects results
#     initial_conditions: Union[dict, monte.ValueMap]
#         Initial conditions for the series. For example:
#
#         .. code-block:: Python
#
#             initial_conditions = {
#                 "temperature": 300.0,
#                 "param_chem_pot": [-4.0],
#             }
#
#     conditions_increment: dict
#         Incremental conditions for the series. For example:
#
#         .. code-block:: Python
#
#             conditions_increment = {
#                 "temperature": 0.0,
#                 "param_chem_pot": [0.5],
#             }
#
#     n_states: int
#         Number of states in the series, including the initial state.
#     output_dir: Union[str, pathlib.Path]
#         Directory to save output files.
#     dependent_runs: bool
#         If True, the initial configuration is set using the final configuration of the
#         previous run.
#     occ_location: Optional[OccLocation] = None
#         Current occupant location list. If provided, the user is
#         responsible for ensuring it is up-to-date with the current
#         occupation of `state`. It is used and updated during the run.
#         If None, an occupant location list is generated for the run.
#
#     Returns
#     -------
#     run_manager: RunManager
#         The input `run_manager` with collected results.
#     """
#     output_params = RunDataOutputParams(
#         do_save_all_initial_states=True,
#         do_save_all_final_states=True,
#         write_initial_states=True,
#         write_final_states=True,
#         output_dir=pathlib.Path(output_dir),
#     )
#     config_generator = FixedConfigGenerator(
#         supercell=supercell,
#         motif=motif,
#     )
#
#     state_generator = IncrementalConditionsStateGenerator(
#         output_params=output_params,
#         initial_conditions=initial_conditions,
#         conditions_increment=conditions_increment,
#         n_states=n_states,
#         config_generator=config_generator,
#         dependent_runs=dependent_runs,
#         modifiers=[],
#     )
#
#     log = MethodLog()
#     log.reset_to_stdout()
#     log.section("Begin: Monte Carlo calculation series")
#
#     log.print("Checking for completed runs...\n")
#     state_generator.read_completed_runs()
#     log.print(f"Found {state_generator.n_completed_runs}\n\n")
#
#     while not state_generator.is_complete:
#         run_manager.run_index = state_generator.n_completed_runs + 1
#
#         log.print("Generating next state...\n")
#         state = state_generator.next_state
#         log.print(xtal.pretty_json(state.conditions.to_dict()))
#         log.print("Done\n")
#
#         run_data = RunData(
#             initial_state=state,
#         )
#
#         log.print(f"Performing Run {run_manager.run_index}...\n")
#         calculator.run(
#             state=state,
#             run_manager=run_manager,
#         )
#         log.print(f"Run {run_manager.run_index} Done\n\n")
#
#         run_data.final_state = state
#         state_generator.append(run_data)
#         state_generator.write_completed_runs()
#     log.print("Monte Carlo calculation series complete\n")
