R"""CASM cluster expansion semi-grand canonical Monte Carlo

Includes:

- the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalConditions`
  class for representing thermodynamic conditions,
- the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalEventGenerator`
  class for proposing events in the semi-grand canonical ensemble (currently only used
  internally),
- the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalPotential`
  class for calculating changes in the semi-grand canonical energy due to the
  proposed events, and
- the :class:`~libcasm.clexmonte.semigrand_canonical.SemiGrandCanonicalCalculator`
  class for sampling microstates in the semi-grand canonical ensemble.

It also makes use of:

- :class:`~libcasm.monte.sampling.SamplingFixture` and
  :class:`~libcasm.monte.sampling.RunManager`, to control sampling, convergence
  checking, and results output.

.. rubric:: Standard sampling functions:

- ``"formation_energy"``: (default=sampled) - Formation energy,
  :math:`E_\mathbb{C}/N_u` (eV), where :math:`E_\mathbb{C}` is the formation
  energy of the microstate :math:`\mathbb{C}`, as calculated by the
  cluster expansion Hamiltonian calculator with key ``"formation_energy"`` obtained from
  :py:attr:`System.clex <libcasm.clexmonte.System.clex>`, and :math:`N_u` is the number
  of unit cells in the supercell.
- ``"mol_composition"``: (default=sampled) - Sample composition per unit cell,
  :math:`\vec{n}_\mathbb{C}`, an `s`-dimensional vector, where `s` is the number of
  species in the system and :math:`n_i = N_i/N_u`, is the number of the `i`-th
  component per unit cell, :math:`N_i` being the total number of the `i`-th component
  (per supercell).
- ``"param_composition"``: (default=sampled) - Sample parametric composition,
  :math:`\vec{x}_\mathbb{C}`, a `k`-dimensional vector, defined as

  .. math::

      \vec{n} = \vec{n}_0 + \pmb{Q} \vec{x},

  with respect to a choice of origin, :math:`\vec{n}_0`, and composition axes,
  :math:`\pmb{Q}=\left[\vec{q}_1, \dots, \vec{q}_k \right]` which span the
  `k`-dimensional space of allowed compositions. The matrix :math:`\pmb{Q}`, has
  dimensions :math:`s \times k`, where `k < s`, due to the fixed number of sites in the
  semi-grand canonical ensemble.
- ``"potential_energy"``: (default=sampled) - Potential energy,
  :math:`\Omega_\mathbb{C}/N_u` (eV), where

  .. math::

      \Omega_\mathbb{C} = E_\mathbb{C} - N_{u}\sum_i\tilde{\mu}_{i} x_{i,\mathbb{C}},

  which determines the microstate probabilities, :math:`P_{\mathbb{C}}`, in the
  semi-grand canonical ensemble according to

  .. math::

      P_{\mathbb{C}} = \
      \frac{e^{-\beta \Omega_\mathbb{C}}}{ \
      \sum_{\mathbb{C}} e^{-\beta \Omega_\mathbb{C}}},

  where :math:`\beta=1/(k_{B}T)`.

- ``"order_parameter_<key>"``: (default=sampled) - Order parameters,
  :math:`\vec{\eta}_\mathbb{C}`, where ``<key>`` is a key in
  :py:attr:`System.dof_space <libcasm.clexmonte.System.dof_space>` specifying
  the :class:`~libcasm.clexulator.DoFSpace` defining the order parameters, as described
  in `Evaluating order parameters`_.
- ``"subspace_order_parameter_<key>"``: (default=included) - Order parameter magnitudes
  by subspace, :math:`\vec{\xi}`, where

  .. math::

      \xi_i = \left({\sum_{j \in \vec{S}_i} \eta_{j,\mathbb{C}}^2}\right)^{1/2},

  :math:`\vec{S}_i` is the set of indices of the order parameter components in the
  `i`-th subspace, and ``<key>`` is a key in
  :py:attr:`System.dof_space <libcasm.clexmonte.System.dof_space>`.
- ``"temperature"``: (default=not sampled) - Temperature, :math:`T` (K). This is a
  constant set by the input conditions.
- ``"param_chem_pot"``: (default=not sampled) - Parametric chemical potential,
  :math:`\vec{\tilde{\mu}}` (eV). This is a constant set by the input conditions.
- ``"formation_energy_corr"``: (default=not sampled) - Formation energy correlation
  functions, :math:`\vec{\Gamma}_{\mathbb{C}}`, which are the mean value of the
  functions in the `i`-th linear function orbit,

  .. math::

      \Gamma_{i,\mathbb{C}} = \frac{1}{m^n_{\alpha}N_u} \
      \sum_{\Phi^m_{\beta} \in \Omega^n_{\alpha}} \Phi^m_{\beta},

  where :math:`m^n_{\alpha}`, is the multiplicity of symmetrically equivalent cluster
  functions, :math:`\Phi^m_{\beta}`, per primitive cell in the orbit
  :math:`\Omega^n_{\alpha}` which corresponds to `i`. The indices `\alpha` and `\beta`
  indicate clusters of sites, and the indices `m` and `n` indicate functions on the
  respective clusters.


.. rubric:: Standard JSON sampling functions (not sampled by default):

- ``"config"``: (default=not sampled) - Store a JSON representation of the
  configuration, using
  :func:`MonteCarloConfiguration.to_dict \
  <libcasm.clexmonte.MonteCarloConfiguration.to_dict>`.


.. rubric:: Standard analysis functions:

- ``"heat_capacity"``: (default=included) - Heat capacity,

  .. math::

        c_{P,\vec{\tilde{\mu}},N_{u}}=\
        \frac{\langle \Omega_i \Omega_j \rangle - \
        \langle \Omega_i \rangle \langle \Omega_j \rangle}{N_{u}k_{B}T^2},

  :math:`k_{B}` being the Boltzmann constant, ``8.6173303E-05`` (eV/K).

- ``"param_susc"``: (default=included) - Parametric chemical susceptibility,

  .. math::

      \chi_{i,j}=\
      \frac{\langle x_i x_j \rangle - \
      \langle x_i \rangle \langle x_j \rangle}{k_{B}T}N_{u}

- ``"mol_susc"``: (default=included) - Mol chemical susceptibility,

  .. math::

      \chi^n_{i,j}=\
      \frac{\langle n_i n_j \rangle - \
      \langle n_i \rangle \langle n_j \rangle}{k_{B}T}N_{u}

- ``"param_thermochem_susc"``: (default=included) - Parametric thermochemical
  susceptibility,

  .. math::

      \chi_{i,\Omega}=\
      \frac{\langle x_i \Omega \rangle - \
      \langle x_i \rangle \langle \Omega \rangle}{k_{B}T}

- ``"mol_thermochem_susc"``: (default=included) - Mol thermochemical susceptibility,

  .. math::

      \chi^n_{i,\Omega}=\
      \frac{\langle n_i \Omega \rangle - \
      \langle n_i \rangle \langle \Omega \rangle}{k_{B}T}


.. _Evaluating order parameters: https://prisms-center.github.io/CASMcode_pydocs/libcasm/clexulator/2.0/usage/order_parameters.html#evaluating-order-parameters

"""

from ._clexmonte_semigrand_canonical import (
    SemiGrandCanonicalConditions,
    SemiGrandCanonicalPotential,
)
from ._SemiGrandCanonicalCalculator import (
    SemiGrandCanonicalCalculator,
    make_initial_state,
)
