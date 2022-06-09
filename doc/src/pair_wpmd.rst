.. index:: pair_style wpmd/cut

pair_style wpmd/cut command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style wpmd/cut Rc

* Rc = global cutoff, -1 means cutoff of half the shortest box length

Examples
""""""""

.. code-block:: LAMMPS

   pair_style wpmd/cut -1
   pair_style wpmd/cut 40.0

Description
"""""""""""

This pair style represents the basic implementation of the Wave Packet Molecular Dynamics (WPMD) :ref:`(Klakow, 1994) <Klakow1994>` for studying nonideal (strongly coupled) systems of charged particles such as the nonideal plasma and warm dense matter. This method is an extension of the classical molecular dynamics (MD) of electrons and ions, where the ions (nuclei) are treated as classical point-like particles and the electrons are represented as normalized Gaussian wavepackets with dynamical width (size). It allows for studying equilibrium states and non-equilibrium processes beyond the Born-Oppenheimer approach due to the explicit dynamics of electrons. At the moment the method is verified for hydrogen and helium plasmas although it is expected to be applicable for heavier atoms as well.

In this method, the single-electron wave function is parametrized by a set of eight time-dependent scalars: the wavepacket position :math:`\mathbf{r}` (3d vector), the wavepacket width :math:`s` (scalar) and their conjugate momenta :math:`\mathbf{p}` (3d vector), :math:`p_s` (scalar):

.. math::

  \varphi(\mathbf{x}) = \left( \frac{3}{2\pi s^2} \right)^{3/4}
  \exp \left\{
  - \left(\frac{3}{4s^2} - \frac{\mathrm{i}{p_s}}{2\hbar s} \right)
  (\mathbf{x}-\mathbf{r})^2 + \frac{\mathrm{i}}{\hbar}{\mathbf{p}}
  \cdot (\mathbf{x}-\mathbf{r})
  \right\}.

The electron Force Field (eFF) (see :doc:`pair_style eff/cut <pair_eff>` was the first pair style of such kind implemented in LAMMPS  being in fact an extension of the original WPMD algorithm where the spin-dependent Pauli potential is added (see below). The definition of the wavepacket width (size) :math:`s` in eFF differs from this pair style by a factor of :math:`\sqrt{3}`.

Within the Hartree approximation the many-electron wave function is given as

.. math::

  \Psi(\{\mathbf{x}_k\}) = \prod_{k=1}^{N_\mathrm{e}}\varphi(\mathbf{x}_k),

where :math:`N_\mathrm{e}` is the number of electrons.

The WPMD model can be used to perform either Mote-Carlo sampling or MD simulations. In both cases the total energy of the system is given by the Hamiltonian

.. math::

  H_\mathrm{wpmd}
  = \left\langle \Psi \right| \hat{H}_\mathrm{wpmd} \left| \Psi \right\rangle
  = K_\mathrm{i} + K_\mathrm{e} + K'_\mathrm{e} + U_\mathrm{ii}
  + U_\mathrm{ei} + U_\mathrm{ee} + U_\mathrm{ext},

where :math:`K_\mathrm{i}` and :math:`K_\mathrm{e} + K'_\mathrm{e}` are the kinetic energies electrons and ions, :math:`U_\mathrm{ii}`, :math:`U_\mathrm{ei}`, :math:`U_\mathrm{ee}` are the potential energies of ion-ion, electron-ion and electron-electron interactions respectively:

.. math::

  & K_\mathrm{i} = \sum\limits_{k=1}^{N_\mathrm{i}}   \frac{\mathbf{p_\mathrm{i}}_k^2}{2m_\mathrm{i}}, \\
  & K_\mathrm{e} = \sum_{k=1}^{N_\mathrm{e}} \frac{\mathbf{p}_k^2}{2m_\mathrm{e}}, \\
  & K'_\mathrm{e} = \sum_{k=1}^{N_\mathrm{e}} \left( \frac{9\hbar^2}{8m_\mathrm{e}s^2_k} + \frac{{p_s}_k^2}{2m_\mathrm{e}} \right), \\
  & U_\mathrm{ii} = \sum\limits_{k<l}^{N_\mathrm{i},\,N_\mathrm{i}}\! \frac{Z^2 e^2}{\left| \mathbf{R}_k - \mathbf{R}_l \right|}, \\
  & U_\mathrm{ei} = - \sum_{k,l}^{N_\mathrm{e},\,N_\mathrm{i}}\! \frac{Z e^2}{|\mathbf{r}_k-\mathbf{R}_l|}\, \mathrm{erf} \!\Bigg(\frac{\sqrt{3}|\mathbf{r}_k - \mathbf{R}_l|}{\sqrt{2}s_k}\Bigg), \\
  & U_\mathrm{ee} = \sum_{k<l}^{N_\mathrm{e},\,N_\mathrm{e}}\! \frac{e^2}{|\mathbf{r}_k - \mathbf{r}_l|}\, \mathrm{erf} \!\Bigg(\frac{\sqrt{3}|\mathbf{r}_k - \mathbf{r}_l|}{\sqrt{2}(s^2_k+s^2_l)^{1/2}}\Bigg), \\

:math:`N_\mathrm{i}` is the number of ions, :math:`m_\mathrm{i}` and :math:`Ze` are the mass and charge of the ions, :math:`\mathbf{R}_k` and :math:`{p_\mathrm{i}}_k` are the position and momentum of ions, :math:`m_\mathrm{e}` and :math:`e` are the electron mass and charge, :math:`U_\mathrm{ext}` is an external potential, e.g.\ the wall boundary. Note that although :math:`K'_\mathrm{e}` is the kinetic energy, in the log and dump files, it is assigned to the potential energy in order to keep the definition of kinetic energy of electron :math:`K_\mathrm{e}` similar to the classical system.

For MD simulations the equations of motion follow from the time-dependent Schrodinger equation. In the case of Hartree approximation, they correspond to the Hamiltonian equations where the width of each electron represents an additional degree of freedom. The Monte-Carlo algorithm is also similar to those of the classical system but involves also the variation of the wavepacket widths.

The time-dependent dynamics of wavepackets is implemented using :doc:`fix nve/wpmd <fix_nve_wpmd>` for the energy conservative system and  \href{fix_nh_wpmd.html}{fix nvt/wpmd, npt/wpmd, nph/wpmd} for the canonical, isothermal-isobaric, and isenthalpic ensembles. The Monte-Carlo sampling is given by :doc:`fix mc/wpmd <fix_mc_wpmd>`. It can be used for all WPMD modifications listed below.

The pair style has only one parameter: Rc is the cutoff radius for the Coulomb interaction.

The pair style is designed to be used with :doc:`atom_style wavepacket <atom_style>` definitions to handle the description of systems with interacting ions and explicit electrons.

There are a few modifications of the original WPMD algorithm to account for the antisymmetry in the many-electron wave function which is related to the exchange-correlation effects, electron degeneracy, and Pauli blocking. In all these methods the spin projections (up or down) are constantly attributed to all electrons. Below is the list of such modifications included in LAMMPS:

* Electron Force Field (eFF), see :doc:`pair_style eff/cut <pair_eff>`. The antisymmetry is implemented via the spin-dependent semi-empirical Pauli potential. The method is almost as fast as the original WPMD.

* Antisymmetrized Wave Packet Molecular Dynamics (AWPMD), see :doc:`pair_style awpmd/cut <pair_awpmd>`. In AWPMD, the many-body wave function for electrons with the same spin projection is antisymmetrized as defined by the unrestricted Hartree-Fock approximation. This pair style also supports the representation of a single electron by multiple Gaussian wavepackets which improves the accuracy of electron-ion bound states (atoms and molecules). Due to computations of the norm-matrix and the modified equations of motion, this method is slower than the original WPMD or eFF.

* The joint method of the Wave Packet Molecular Dynamics and the Density Functional Theory (WPMD-DFT), see :doc:`pair_style wpmd/dft/cut <pair_wpmd_dft>`. In this method, the exchange-correlation effects are included by an additional energy term computed as a functional of the local electron density following the idea of DFT within the LSDA approximation. The electron density is obtained from the wavepacket positions and widths. This method is slower than eFF but faster than AWPMD. It has a GPU-accelerated version.

The original WPMD method has a known problem of unlimited broadening of the wavepackets for weakly bound electrons. Therefore the periodic boundaries are appropriate only for very high electron density. This problem can be solved either by manual limiting of the wavepacket width (eFF) or by using a 3-dimensional confining potential (wall potential) for the whole system which naturally constrains both the wavepacket positions and widths (see :doc:`fix wall/wpmd <fix_wall_wpmd>`).


Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for this pair style. For minimization of the system's energy with `wpmd/cut` pair
should use a :doc:`fix mc/wpmd <fix_mc_wpmd>`.

This pair style writes its information to :doc:`binary restart files <restart>`,
so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""
This pair is part of the AWPMD package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

This pair is work only with *real* and *electron* units due to energy conversion units.

This pair required :doc:`wavepacket <atom_style>` or :doc:`electron <atom_style>` atom style.

The time evolution can be calculated by one of :doc:`fix nvt/wpmd <fix_nh_wpmd>`, :doc:`fix nph/wpmd <fix_nh_wpmd>`, :doc:`fix npt/wpmd <fix_nh_wpmd>` or :doc:`fix nve/wpmd <fix_nve_wpmd>`. The Monte-Carlo sampling can be performed with :doc:`fix mc/wpmd <fix_mc_wpmd>`.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

:doc:`pair_style eff/cut <pair_eff>`

:doc:`pair_style awpmd/cut <pair_awpmd>`

:doc:`pair_style wpmd/dft/cut <pair_wpmd_dft>`

----------

.. _Klakow1994:

**(Klakow, 1994)** D. Klakow, C. Toepffer, and P.-G. Reinhard, Semiclassical molecular dynamics for strongly coupled coulomb systems, J. Chem. Phys., 101:10766 (1994).
