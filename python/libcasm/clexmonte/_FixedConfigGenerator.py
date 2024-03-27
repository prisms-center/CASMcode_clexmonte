import copy

from ._clexmonte_state import (
    MonteCarloConfiguration,
)
from ._RunData import (
    RunData,
)


class FixedConfigGenerator:
    """A `ConfigGenerator` for state generation - always returns the same configuration

    Notes
    -----

    - Returns the same configuration no matter what the current conditions and
      completed runs are.

    """

    def __init__(
        self,
        configuration: MonteCarloConfiguration,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        configuration: MonteCarloConfiguration
            The MonteCarloConfiguration to generate.
        """
        self._configuration = copy.copy(configuration)

    def __call__(
        self,
        completed_runs: list[RunData],
    ) -> MonteCarloConfiguration:
        """Always returns the same configuration

        Parameters
        ----------
        completed_runs: list[RunData]
            Ignored, but included for compatibility.

        Returns
        -------
        configuration: MonteCarloConfiguration
            The MonteCarloConfiguration provided to the constructor.
        """
        return copy.copy(self._configuration)
