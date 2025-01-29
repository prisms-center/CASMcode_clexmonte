"""Helpers for iterating over sites in a Monte Carlo simulation

The site iterators in this module yield :class:`~libcasm.clexmonte.SiteInfo` that are
updated at each step to give the position, coordinates, etc. for a particular set of
sites in a Monte Carlo simulation. This can be used to help implement custom sampling
and analysis functions.

.. rubric:: Example: Sampling the local composition around vacancies

.. code-block:: python

    import libcasm.clexmonte as clexmonte
    import libcasm.monte.sampling as sampling

    # calculator: clexmonte.MonteCalculator
    # sampling_fixture: clexmonte.SamplingFixture

    va_site_iterator = clexmonte.ChemicalTypeSiteIterator(
        chemical_name=["Va"],
        calculator=calculator,
    )

    def my_custom_function():
        configuration = calculator.state_data.state.configuration

        for site in va_site_iterator:
            # print("site.l:", site.l, "occ: ", configuration.occ(site.l))

        ... do something ...

        return [component_value1, ...]

    sampling_f = sampling.StateSamplingFunction(
        name="my_custom_function",
        description="desc",
        shape=[],
        function=my_custom_function,
        component_names=["component_name1", ...],
    )

    sampling_fixture.add_sampling_function(sampling_f)
    sampling_fixture.sample("my_custom_function")

"""

import typing

import libcasm.clexmonte as clexmonte
import libcasm.monte.events as monte_events


def chemical_name_to_candidate_indices(
    chemical_name: typing.Union[str, list[str]],
    occ_location: monte_events.OccLocation,
):
    """Return a list of all OccCandidate indices for a given chemical species

    Parameters
    ----------
    chemical_name : Union[str, list[str]]
        The chemical name (as given in
        :func:`Occupant.name <libcasm.xtal.Occupant.name>`) of one or more
        occupants.
    occ_location : monte_events.OccLocation
        The occupant location list.

    Returns
    -------
    candidate_indices : list[int]
        Indices into `occ_location.candidate_list()` corresponding to the specified
        chemical species on all asymmetric unit sites.
    """
    if isinstance(chemical_name, str):
        chemical_name = [chemical_name]

    convert = occ_location.convert()
    occ_candidate_list = occ_location.candidate_list()
    candidate_indices = []

    for species_index in range(convert.species_size()):
        occupant = convert.species_index_to_occupant(species_index)

        if occupant.name() not in chemical_name:
            continue

        for candidate_index, cand in enumerate(occ_candidate_list):
            if cand.species_index == species_index:
                candidate_indices.append(candidate_index)
    return candidate_indices


def occupant_name_to_candidate_indices(
    occupant_name: typing.Union[str, list[str]],
    occ_location: monte_events.OccLocation,
):
    """Return a list of all OccCandidate indices for a given occupant

    Parameters
    ----------
    occupant_name : Union[str, list[str]]
        The occupant DoF names (as given in
        :func:`Prim.occ_dof <libcasm.xtal.Prim.occ_dof>`) of one or more occupants.
    occ_location : monte_events.OccLocation
        The occupant location list.

    Returns
    -------
    candidate_indices : list[int]
        Indices into `occ_location.candidate_list()` corresponding to the specified
        occupant(s) on all asymmetric unit sites.
    """
    if isinstance(occupant_name, str):
        occupant_name = [occupant_name]

    convert = occ_location.convert()
    occ_candidate_list = occ_location.candidate_list()
    candidate_indices = []

    for _name in occupant_name:
        species_index = convert.species_name_to_index(species_name=_name)

        for candidate_index, cand in enumerate(occ_candidate_list):
            if cand.species_index == species_index:
                candidate_indices.append(candidate_index)
    return candidate_indices


class SiteInfo:
    """Convenience class for obtaining information about a site and its occupant"""

    def __init__(
        self,
        occ_location: monte_events.OccLocation,
        l: typing.Optional[int] = None,
        mol_id: typing.Optional[int] = None,
    ):
        """
        .. rubric:: Constructor

        Parameters
        ----------
        occ_location: libcasm.monte.events.OccLocation
            The occupant location list.
        l: Optional[int] = None
            The linear site index of the current site, :math:`l`. One and only one of
            `l` or `mol_id` must be provided.
        mol_id: Optional[int] = None
            Index of the :class:`~libcasm.monte.events.Mol` on the current site in the
            occupant location `mol` list. One and only one of `l` or `mol_id` must be
            provided.
        """

        self.occ_location = occ_location
        """libcasm.monte.events.OccLocation: The occupant location list."""

        self.convert = occ_location.convert()
        """libcasm.monte.events.Conversions: Data structure used for index 
        conversions.
        """

        self.l = None
        """int: The linear site index of the current site, :math:`l`.
        
        This is equivalent to py:attr:`SiteInfo.linear_site_index`.
        """

        self.mol_id = None
        """int: Index of the :class:`~libcasm.monte.events.Mol` on the current site in 
        the occupant location `mol` list.
        """

        if l is not None:
            self.set_l(l=l)
        elif mol_id is not None:
            self.set_mol_id(mol_id=mol_id)
        else:
            raise ValueError(
                "Error in SiteInfo constructor: One and only one of `l` or `mol_id` "
                "must be provided."
            )

    def set_l(
        self,
        l: int,
    ):
        """Set the linear site index, :math:`l`, for the current site, and use it to
        set `mol_id`.

        Parameters
        ----------
        l : int
            The linear site index of the current site, :math:`l`.

        """
        self.l = l
        self.mol_id = self.occ_location.linear_site_index_to_mol_id(l)

    def set_mol_id(self, mol_id: int):
        """Set the `mol_id` for the current site, and use it to set the linear
        site index, :math:`l`).

        Parameters
        ----------
        mol_id : int
            Index of the :class:`~libcasm.monte.events.Mol` in the
            occupant location `mol` list.

        """
        self.mol_id = mol_id
        self.l = self.occ_location.mol(mol_id).linear_site_index

    def __copy__(self):
        return SiteInfo(
            mol_id=self.mol_id,
            occ_location=self.occ_location,
        )

    def __deepcopy__(self, memo):
        return SiteInfo(
            mol_id=self.mol_id,
            occ_location=self.occ_location,
        )

    def copy(self):
        return self.__copy__()

    @property
    def linear_site_index(self):
        """int: The linear site index of the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.l`.
        """
        return self.l

    @property
    def asym(self):
        """int: The asymmetric unit index, :math:`a`, of the current site, :math:`l`."""
        return self.convert.l_to_asym(self.l)

    @property
    def b(self):
        """int: The sublattice / basis index of the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.sublattice`.
        """
        return self.convert.l_to_b(self.l)

    @property
    def sublattice(self):
        """int: The sublattice / basis index of the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.b`.
        """
        return self.convert.l_to_b(self.l)

    @property
    def basis_cart(self):
        """numpy.ndarray[numpy.float64[3, 1]]: The Cartesian coordinate,
        :math:`r_{cart}`, in the primitive unit cell, of the sublattice that the
        current site, :math:`l`, belongs to."""
        return self.convert.l_to_basis_cart(self.l)

    @property
    def basis_frac(self):
        """numpy.ndarray[numpy.float64[3, 1]]: The fractional coordinate,
        :math:`r_{frac}`, in the primitive unit cell, of the sublattice that the
        current site, :math:`l`, belongs to."""
        return self.convert.l_to_basis_frac(self.l)

    @property
    def bijk(self):
        """libcasm.xtal.IntegralSiteCoordinate: The integral site coordinates,
        :math:`(b,i,j,k)` for the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.integral_site_coordinate`.
        """
        return self.convert.l_to_bijk(self.l)

    @property
    def integral_site_coordinate(self):
        """libcasm.xtal.IntegralSiteCoordinate: The integral site coordinates,
        :math:`(b,i,j,k)` for the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.bijk`.
        """
        return self.convert.l_to_bijk(self.l)

    @property
    def ijk(self):
        """numpy.ndarray[numpy.int64[3, 1]]: The unit cell indices, :math:`(i,j,k)`,
        for the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.unitcell_indices`.
        """
        return self.convert.l_to_ijk(self.l)

    @property
    def unitcell_indices(self):
        """numpy.ndarray[numpy.int64[3, 1]]: The unit cell indices, :math:`(i,j,k)`,
        for the current site, :math:`l`.

        This is equivalent to py:attr:`SiteInfo.ijk`.
        """
        return self.convert.l_to_ijk(self.l)

    @property
    def linear_unitcell_index(self):
        """int: The linear unitcell index"""

    @property
    def coordinate_cart(self):
        """numpy.ndarray[numpy.float64[3, 1]]: The Cartesian coordinate,
        :math:`r_{cart}`, of the current site,: math:`l`.
        """
        return self.convert.l_to_cart(self.l)

    @property
    def coordinate_frac(self):
        """numpy.ndarray[numpy.float64[3, 1]]: The fractional coordinate,
        :math:`r_{frac}`, relative to the :class:`~libcasm.xtal.Prim` lattice vectors,
        :math:`P`, from the linear site index,: math:`l`.
        """
        return self.convert.l_to_frac(self.l)

    @property
    def unitl(self):
        """int: The non-primitive unit cell sublattice index, :math:`l'`, of the
        current site, :math:`l`.

        Linear site index in a non-primitive unit cell. When a non-primitive unit cell
        is used to construct a supercell and set the appropriate symmetry for a problem,
        conversions between :math:`l`, :math:`b`, and :math:`l'` may all be useful.

        This is primarily a placeholder for possible future use cases in which the
        asymmetric unit is not contained in the primitive unit cell. In the typical
        case, where the unit cell is the primitive cell, then :math:`l' == b`.
        """
        return self.convert.l_to_unitl(self.l)

    @property
    def mol(self):
        """libcasm.monte.events.Mol: The :class:`~libcasm.monte.events.Mol` object for
        the current site, :math:`l`.
        """
        return self.occ_location.mol(self.mol_id)

    @property
    def occupant(self):
        """libcasm.xtal.Occupant: The current occupant at the current site,
        :math:`l`."""
        return self.convert.species_index_to_occupant(self.mol.species_index)

    @property
    def chemical_name(self):
        """str: The chemical name (as given in
        :func:`Occupant.name <libcasm.xtal.Occupant.name>`) of the occupant at
        the current site, :math:`l`.
        """
        return self.occupant.name()

    @property
    def occupant_name(self):
        """str: The occupant DoF name (as given in
        :func:`Prim.occupants <libcasm.xtal.Prim.occupants>`) of the occupant at
        the current site, :math:`l`.
        """
        return self.convert.species_index_to_name(self.mol.species_index)

    def _atoms(self):
        """Return a list of atoms at the current site, :math:`l`."""
        atoms = []
        for atom_id in self.mol.component_id:
            atoms.append(self.occ_location.atom(atom_id))
        return atoms

    @property
    def atom_id(self):
        """int: (Applies to KMC simulations with atomic occupants) The atom id at the
        current site, :math:`l`.

        The `atom_id` is the index into the current list of atoms,
        :func:`OccLoation.atom <libcasm.monte.events.OccLocation.atom>`.
        """
        return self.mol.component_id[0]

    @property
    def unique_atom_id(self):
        """int: (Applies to KMC simulations with atomic occupants) The unique atom id
        for the atom at the current site, :math:`l`.

        The `unique_atom_id` differs from the `atom_id` for KMC simulations which allow
        deposition / dissolution type events and is the key used for the
        :func:`OccLoation.atom_info_initial
        <libcasm.monte.events.OccLocation.atom_info_initial>`
        and
        :func:`OccLoation.atom_info_final
        <libcasm.monte.events.OccLocation.atom_info_final>`
        data.
        """
        return self.occ_location.unique_atom_id[self.atom_id]

    @property
    def atom(self):
        """libcasm.monte.events.Atom: (Applies to KMC simulations with atomic
        occupants) The atom at the current site, :math:`l`."""
        return self.occ_location.atom(self.mol.component_id[0])

    @property
    def atom_n_jumps(self):
        """int: (Applies to KMC simulations with atomic
        occupants) The atom at the current site, :math:`l`."""
        return self.atom.n_jumps

    @property
    def atom_coordinate_frac(self):
        """numpy.ndarray[numpy.float64[3, 1]]: (Applies to KMC simulations with atomic
        occupants) Current position, in Cartesian coordinates, as if no periodic
        boundary, for the occupant of the current site, :math:`l`.

        Note that this is not a displacement vector, but the absolute position of the
        occupant in the current configuration, if no periodic boundary conditions were
        applied.
        """
        return self.atom.translation + self.basis_frac

    @property
    def atom_coordinate_cart(self):
        """numpy.ndarray[numpy.float64[3, 1]]: (Applies to KMC simulations with atomic
        occupants) Current position, in Cartesian coordinates, as if no periodic
        boundary, for the occupant of the current site, :math:`l`.

        Note that this is not a displacement vector, but the absolute position of the
        occupant in the current configuration, if no periodic boundary conditions were
        applied.
        """
        return self.convert.lat_column_mat() @ self.atom_coordinate_frac

    @property
    def component_id(self):
        """LongVector: (Applies to KMC simulations with multi-atom occupants) A list of
        atom id for the occupant at the current site, :math:`l`.

        The `atom_id` is the index into the current list of atoms,
        :func:`OccLoation.atom <libcasm.monte.events.OccLocation.atom>`."""
        return self.mol.component_id

    @property
    def component_atoms(self):
        """list[libcasm.monte.events.Atom]: (Applies to KMC simulations with multi-atom
        occupants) The component atoms of the occupant at the current site, :math:`l`.
        """
        return [self.occ_location.atom(atom_id) for atom_id in self.mol.component_id]


class CompleteSiteIterator:
    """Yield :class:`~libcasm.clexmonte.SiteInfo` for all sites in the system"""

    def __init__(
        self,
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """
        ..rubric:: Constructor

        Parameters
        ----------
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list. One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested. One and only one of
            `occ_location` and `calculator` must not be `None`.
        """
        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in CompleteSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for all sites in the system"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        self._site_info = SiteInfo(
            l=0,
            occ_location=self.occ_location,
        )

        convert = self.occ_location.convert()

        for l in range(convert.l_size()):
            self._site_info.set_l(l=l)
            yield self._site_info


class SublatticeSiteIterator:
    """Yield :class:`~libcasm.clexmonte.SiteInfo` for 1 or more sublattices"""

    def __init__(
        self,
        sublattice: typing.Union[int, list[int]],
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """
        ..rubric:: Constructor

        Parameters
        ----------
        sublattice: Union[int, list[int]]
            The sublattice index of one or more sublattices.
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list. One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested. One and only one of
            `occ_location` and `calculator` must not be `None`.
        """
        if isinstance(sublattice, int):
            sublattice = [sublattice]

        self.sublattice = sorted(sublattice)
        """list[int]: The sublattice index of one or more sublattices."""

        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in SublatticeSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for the requested sites"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        self._site_info = SiteInfo(
            l=0,
            occ_location=self.occ_location,
        )

        convert = self.occ_location.convert()
        vol = convert.unitcell_index_converter().total_unitcells()

        for b in self.sublattice:
            for l in range(vol * b, vol * (b + 1)):
                self._site_info.set_l(l=l)
                yield self._site_info


class AsymmetricUnitSiteIterator:
    """Yield :class:`~libcasm.clexmonte.SiteInfo` for the sites in one or more
    asymmetric units"""

    def __init__(
        self,
        asym: typing.Union[int, list[int]],
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """
        ..rubric:: Constructor

        Parameters
        ----------
        asym: Union[int, list[int]]
            One or more asymmetric unit indices.
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list. One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested. One and only one of
            `occ_location` and `calculator` must not be `None`.
        """
        if isinstance(asym, int):
            asym = [asym]

        self.asym = sorted(asym)
        """list[int]: Asymmetric unit indices."""

        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in AsymmetricUnitSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for the requested sites"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        convert = self.occ_location.convert()

        for a in self.asym:
            sites_iterator = SublatticeSiteIterator(
                sublattice=convert.asym_to_b(a),
                occ_location=self.occ_location,
            )
            for site_info in sites_iterator:
                yield site_info


class OccCandidateSiteIterator:
    """Yield class: `~libcasm.clexmonte.SiteInfo` for 1 or more OccCandidate"""

    def __init__(
        self,
        candidate_indices: typing.Union[int, list[int]],
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """
        ..rubric:: Constructor

        Parameters
        ----------
        candidate_indices: Union[int, list[int]]
            One or more asymmetric unit indices.
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list. One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested. One and only one of
            `occ_location` and `calculator` must not be `None`.
        """
        if isinstance(candidate_indices, int):
            candidate_indices = list[candidate_indices]

        self.candidate_indices = candidate_indices
        """list[int]: Indices into `occ_location.candidate_list()` corresponding to the 
        type of occupant(s) and site(s) requested."""

        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in SublatticeSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for the requested sites"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        self._site_info = None

        # `candidate_index`: index into list of OccCandidate corresponding to
        # the (species index, asym unit index)
        for candidate_index in self.candidate_indices:
            # Current size of the list occupants corresponding to
            # `candidate_index`-th OccCandidate
            n = self.occ_location.cand_size_by_candidate_index(candidate_index)

            # `location_index`: index into list occupants corresponding to
            # the (species, asym unit) of the `candidate_index`-th OccCandidate
            for location_index in range(n):
                # mol_id: Index into `OccLocation.mol` list
                mol_id = self.occ_location.mol_id_by_candidate_index(
                    candidate_index=candidate_index,
                    location_index=location_index,
                )

                if self._site_info is None:
                    self._site_info = SiteInfo(
                        mol_id=mol_id,
                        occ_location=self.occ_location,
                    )
                else:
                    self._site_info.set_mol_id(mol_id=mol_id)

                yield self._site_info


class ChemicalTypeSiteIterator:
    """Yield :class:`~libcasm.clexmonte.SiteInfo` for 1 or more chemical species"""

    def __init__(
        self,
        chemical_name: typing.Union[str, list[str]],
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """
        ..rubric:: Constructor

        Parameters
        ----------
        chemical_name: Union[str, list[str]]
            The chemical name( as given in
            :func:`Occupant.name <libcasm.xtal.Occupant.name>`) of one or more
            occupants.
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list. One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested. One and only one of
            `occ_location` and `calculator` must not be `None`.
        """

        if isinstance(chemical_name, str):
            chemical_name = [chemical_name]

        self.chemical_name = chemical_name
        """list[str]: The chemical name (as given in: 
        func:`Occupant.name <libcasm.xtal.Occupant.name>`) of one or more occupants.
        """

        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in ChemicalTypeSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for the requested occupant"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        sites_iterator = OccCandidateSiteIterator(
            candidate_indices=chemical_name_to_candidate_indices(
                chemical_name=self.chemical_name,
                occ_location=self.occ_location,
            ),
            occ_location=self.occ_location,
        )

        for site_info in sites_iterator:
            yield site_info


class OccupantTypeSiteIterator:
    """Yield :class:`~libcasm.clexmonte.SiteInfo` for 1 or more occupants

    This yields a :class:`~libcasm.clexmonte.SiteInfo` that is updated at each step
    to give the position, coordinates, etc. for the requested occupant types, based on
    the current occupant location list.
    """

    def __init__(
        self,
        occupant_name: typing.Union[str, list[str]],
        occ_location: typing.Optional[monte_events.OccLocation] = None,
        calculator: typing.Optional[clexmonte.MonteCalculator] = None,
    ):
        """

        ..rubric:: Constructor

        Parameters
        ----------
        occupant_name: Union[str, list[str]]
            The occupant DoF names (as given in
            :func:`Prim.occupants <libcasm.xtal.Prim.occupants>`) of one or more
            occupants.
        occ_location: Optional[monte_events.OccLocation] = None
            The occupant location list.One and only one of `occ_location` and
            `calculator` must not be `None`.
        calculator: Optional[clexmonte.MonteCalculator] = None
            A MonteCalculator which will be used to obtain an occupant location list
            from `calculator.state_data` when requested.One and only one of
            `occ_location` and `calculator` must not be `None`.
        """

        if isinstance(occupant_name, str):
            occupant_name = [occupant_name]

        self.occupant_name = occupant_name
        """list[str]: The occupant DoF names (as given in 
        :func:`Prim.occupants <libcasm.xtal.Prim.occupants>`) of one or more occupants.
        """

        if occ_location is None and calculator is None:
            raise ValueError(
                "Error in OccupantTypeSiteIterator: If `occ_location` is None, "
                "`calculator` must be a MonteCalculator"
            )

        self.occ_location = occ_location
        """Optional[libcasm.monte.events.OccLocation]: The occupant location list."""

        self.calculator = calculator
        """Optional[libcasm.clexmonte.MonteCalculator]: The MonteCalculator used to 
        obtain the occupant location list. """

    def __iter__(self):
        """Yield :class:`~libcasm.clexmonte.SiteInfo` for the requested occupants"""

        if self.occ_location is None:
            self.occ_location = self.calculator.state_data.occ_location

        sites_iterator = OccCandidateSiteIterator(
            candidate_indices=occupant_name_to_candidate_indices(
                occupant_name=self.occupant_name,
                occ_location=self.occ_location,
            ),
            occ_location=self.occ_location,
        )

        for site_info in sites_iterator:
            yield site_info
