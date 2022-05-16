.. index:: pair_style awpmd/cut

pair_style awpmd/cut command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style awpmd/cut Rc keyword value ...

* Rc = global cutoff, -1 means cutoff of half the shortest box length
* zero or more keyword/value pairs may be appended
* keyword = *hartree* or *uhf* or *free* or *pbc* or *fix* or *harm* or *ermscale* or *flex_press*

.. parsed-literal::

       *hartree* value = none
       *uhf* value = none
       *free* value = none
       *pbc* value = Plen
         Plen = periodic width of electron = -1 or positive value (distance units)
       *fix* value = Flen
         Flen = fixed width of electron = -1 or positive value (distance units)
       *harm* value = width
         width = harmonic width constraint
       *ermscale* value = factor
         factor = scaling between electron mass and width variable mass
       *flex_press* value = none

Examples
""""""""

.. code-block:: LAMMPS

   pair_style awpmd/cut -1
   pair_style awpmd/cut 40.0 uhf free
   pair_coeff * *
   pair_coeff 2 2 20.0

Description
"""""""""""

This pair style represents Split AWPMD, the most general version of the Wave Packet Molecular Dynamics with antisymmetrization and wave packet splitting. Compare to the original method of Wave Packet Molecular Dynamics (WPMD) (see :doc:`pair_style wpmd/cut <pair_wpmd>`) it is more precise and accounts for the exchange interaction of electrons although it is about 100 times slower.

In the splitting technique :ref:`(Morozov, 2012) <Morozov2012>` an electron is represented by multiple Gaussians, with mixing coefficients playing the role of additional dynamic variables. It significantly improves the accuracy of the wave function representation and allows one to reproduce many quantum effects such as penetration of a potential barrier. The ions are treated as classical particles.

In the general formulation of wave packet molecular dynamics the many electron trial wave function :math:`\Psi` is parameterized by a set of variables :math:`\mathbf{q}`. Then it is substituted into the time dependent Schroedinger equation and the variational principle is used to obtain equations of motion for the dynamic
variables :math:`\mathbf{q}(t)`:

.. math::

	\mathbf{N}\dot{\mathbf{q}}=\frac{\partial H}{\partial \mathbf{q}},\qquad
	N_{ab}=\frac{\partial}{\partial q^*_a}
	\frac{\partial}{\partial q_b}\ln\langle\Psi(\mathbf{q}^*)|\Psi(\mathbf{q})\rangle.
	
Here :math:`\mathbf{N}` is the norm-matrix, relating the generalized velocities and forces.

In the WPMD/MC method the discretized version of the above equations is solved dynamically or by Monte-Carlo sampling. The Hamiltonian of the many electron system interacting with ions reads:

.. math::

	\hat{H} = \hat{K}^\mathrm{e}+\hat{U}^\mathrm{ei}+\hat{U}^\mathrm{ee}+\hat H_\mathrm{ext}
	= - \sum_k\frac{\hbar^2\Delta_k}{2m}
	- \sum_{k,i}\frac{eq_i}{|\hat{\mathbf{x}}_k-\mathbf{R}_i|}
	+ \sum_{k<m}\frac{e^2}{|\hat{\mathbf{x}}_k-\hat{\mathbf{x}}_m|}
	+ \hat H_\mathrm{ext},

where :math:`\hat{K}^\mathrm{e}` is the electron kinetic energy, :math:`\hat{U}^\mathrm{ei}` is the electron-ion Coulomb interaction, :math:`\hat{U}^\mathrm{ee}` is the electron-electron Coulomb repulsion, and :math:`\hat H_\mathrm{ext}` is the external field potential, :math:`\mathbf{R}_i` and :math:`q_i` are the ion coordinates and charges.

For Split WPMD we expand the single electron wave function :math:`\phi_k(\mathbf{x})` in a number of simple WPs :math:`\varphi_{k\alpha}(\mathbf{x})`:

.. math::

	\phi_k(\mathbf{x}) = n_k^{-1/2}\sum_{\alpha=1}^{M_k} c_{k\alpha} \varphi_{k\alpha}(\mathbf{x}),

with :math:`\varphi_{k\alpha}(\mathbf{x})` being normalized Gaussian wave packets (WPs):

.. math::

	\varphi_{k\alpha}(\mathbf{x}) = \left( \frac{3}{2\pi s_{k\alpha}^2} \right)^{3/4}
	\exp \left\{- \left(\frac{3}{4s_{k\alpha}^2} - \frac{\mathrm{i}{p_{sk\alpha}}}{2\hbar s_{k\alpha}} \right)
	(\mathbf{x}-\mathbf{r}_{k\alpha})^2 + \frac{\mathrm{i}}{\hbar}{\mathbf{p}_{k\alpha}} \cdot (\mathbf{x}-\mathbf{r}_{k\alpha})
	\right\},

and

.. math::

	n_k=\sum_{\alpha,\beta} c^*_{k\alpha}
	c^{}_{k\beta}\int\varphi^*_{k\alpha}\varphi^{}_{k\beta}d^3x

is the normalizing factor for :math:`\phi_k`.

The time-dependent complex coefficients :math:`c_{k\alpha}(t)` together with the standard WP parameters :math:`\mathbf{r}_{k\alpha}(t)`, :math:`\mathbf{p}_{k\alpha}(t)`, :math:`s_{k\alpha}(t)`, :math:`p_{s_{k\alpha}}(t)` constitute the set of dynamic variables for :math:`k`-th electron. As seen,  the variational freedom is extended from 8 real parameters per electron in the original WPMD to 10 parameters in AWPMD. These parameters are seen by LAMMPS via :doc:`atom_style wavepacket <atom_style>`.

The total number of generalized dynamic variables (including both coordinates and conjugate momenta)  is controlled by the number of WPs per electron :math:`M_k`, which may be set different for different electrons. When :math:`M_k=1`, the factor :math:`c` becomes redundant and the scheme reduces to the original WPMD. Classical interpretation of WP parameters (e.g. :math:`\mathbf{r}` being the mean electron position) within the Split WPMD model is valid only for non-overlapping single electron WPs and should be considered with care.

Significant advantage of the Gaussian expansion is that the interaction matrix elements are proportional to the corresponding WP overlaps :math:`o_{k\alpha l\beta}=\int\varphi^*_{k\alpha}\varphi^{}_{l\beta}d^3x`:

.. math::

	&\langle\varphi_{k\alpha}(\mathbf{x})|	-\frac{\hbar^2\Delta}{2m}|	\varphi_{l\beta}(\mathbf{x})\rangle
	= o_{k\alpha l\beta}K^\mathrm{e}_{k\alpha l\beta},\qquad\\
	&\langle\varphi_{k\alpha}(\mathbf{x})|
	-\sum_{i}\frac{eq_i}{|\hat{\mathbf{x}}-\mathbf{R}_i|}|
	\varphi_{l\beta}(\mathbf{x})\rangle
	= o_{k\alpha l\beta}U^\mathrm{ei}_{k\alpha l\beta}, \\
	&\langle\varphi_{k\alpha}(\mathbf{x}_1)\varphi_{l\beta}(\mathbf{x}_2)|
	\frac{e^2}{|\hat{\mathbf{x}}_1-\hat{\mathbf{x}}_2|}|
	\varphi_{m\gamma}(\mathbf{x}_1)\varphi_{n\delta}(\mathbf{x}_2)\rangle
	= o_{k\alpha m\gamma}o_{l\beta n\gamma}U^\mathrm{ee}_{k\alpha l\beta m\gamma n\delta}.

The WP overlaps :math:`o_{k\alpha l\beta}` and residual matrix elements :math:`K^\mathrm{e}_{k\alpha l\beta}`, :math:`U^\mathrm{ei}_{k\alpha l\beta}`, :math:`U^\mathrm{ee}_{k\alpha l\beta	m\gamma n\delta}` are easily obtained analytically for the Gaussian WPs as well as their derivatives with respect to the WP parameters (see :ref:`(Valuev, 2015) <Valuev2015>`).

The total many-electron wave function may be constructed by using different quantum approximations accounting for electron spins. The most used for WPMD are the Hartree approximation (trial state is the product of single electron wave functions) :ref:`(Klakow, 1994) <Klakow1994awpmd>` and the antisymmetrized approximation :ref:`(Jakob, 2007) <Jakob2007>`. The latter is equivalent to the unrestricted Hartree-Fock (UHF) approach when trial state is a single determinant of spin orbitals. The spin orbitals are constructed by explicitly associating spin up or spin down state with each of the spatial single electron wave functions :math:`\phi_k`. The total energy for both Hartree and UHF cases in the Split WPMD model reads:

.. math::

      H &= \!\!\!\!\sum_{\small {\rm same~spin}~(k, l)}\!\!
(n_kn_l)^{-\frac{1}{2}}y_{kl}\sum_{\alpha,\beta}(c^*_{k\alpha}o^{}_{k\alpha l\beta}c^{}_{l\beta})
      (K^e_{k\alpha l\beta}+U^\mathrm{ei}_{k\alpha l\beta})  \\
      & {} + \!\!\!\!\sum_{\small {\rm same~spin}~(k, m)}\sum_{\small {\rm same~spin}~(l, n) }\!\!
      (n_kn_ln_mn_n)^{-\frac{1}{2}}y_{mk}y_{nl}\sum_{\alpha,\beta,\gamma,\delta}
      (c^*_{k\alpha}o^{}_{k\alpha m\gamma}c^{}_{m\gamma})(c^*_{l\beta}o^{}_{l\beta n\delta}c^{}_{n\delta})
      U^\mathrm{ee}_{k\alpha l\beta m\gamma n\delta} \\
      & {} - k_\mathrm{exch}\!\!\sum_{\small{\rm same~spin}~(k, l, m, n)}\!\!
      (n_kn_ln_mn_n)^{-\frac{1}{2}}y_{ml}y_{nk}\sum_{\alpha,\beta,\gamma,\delta}
      (c^*_{k\alpha}o_{k\alpha m\gamma}c_{m\gamma})(c^*_{l\beta}o_{l\beta n\delta}c_{n\delta})
      U^\mathrm{ee}_{k\alpha l\beta m\gamma n\delta}.

For the UHF case (AWPMD) :math:`k_\mathrm{exch}=1` and :math:`y_{ij}` are elements of the inverse overlap matrix :math:`\mathbf{Y} = \mathbf{O}^{-1}` (:math:`O_{km}=\sum_{\alpha,\beta}c^*_{k\alpha}o^{}_{k\alpha m\beta}c^{}_{m\beta}` for same spin orbital indices :math:`k` and :math:`m`). In this pair style it corresponds to the keyword *uhf*.

The Hartree case is recovered by setting :math:`k_\mathrm{exch}=0` and :math:`y_{ij}=\delta_{ij}`, thus greatly simplifying the summation. It corresponds to the keyword *hartree*. Note that when the number of WPs per electron is 1, the *hartree* setting should produce the same energy as :doc:`pair_style wpmd/cut <pair_wpmd>`, however the computation speed would be significantly lower. 

Wave packet width restriction is an important part of wave packet simulation. It may be either controlled by the internal pair style settings or by applying a constraint potential walls to the whole system (see :doc:`fix wall/wpmd <fix_wall_wpmd>`).

The *free*, *pbc*, *fix* and *harm* keywords specify the internal constraints on the electron wave width.  

If the *free* keyword is specified, then there is no width constraint. This setting is default and is useful for the ground state or low temperature computations. Also *free* setting should be used in conjunction with the constraint potential walls provided by :doc:`fix wall/wpmd <fix_wall_wpmd>` (recommended).  

If the *fix* keyword is used and *Flen* is specified as -1, then wave packets have a constant widths that are read from the data file or kept at their initial values after construction with :doc:`create atoms <create_atoms>`.

The simplest dynamic restriction are periodical boundary conditions for the WP widths. In this mode the widths can grow only up to some maximum value *Plen*, any larger :math:`s` values are treated as 2Plen-s and the width momentum is reversed. If the *pbc* keyword is used and *Plen* is specified as -1, then the maximum width is half the shortest box length.  If *Plen* is a positive value, then the value is the maximum width. Note that the periodic boundary conditions do not solve the broadening problem completely and the simulation results usually depend on *Plen*. 

A more elaborate solution proposed in :ref:`(Zwicknagel, 2006) <Zwicknagel2006>` is to introduce an additional harmonic term :math:`\Delta H = (9\hbar^2s_k^2)/(8ms_0^4)` to the Hamiltonian which prevents WP from spreading. It corresponds to the *harm* keyword. The free parameter :math:`s_0` stands for the mean value of width :math:`s` if there were no Coulomb interaction. In :ref:`(Zwicknagel, 2006) <Zwicknagel2006>` it was taken to be :math:`s_0 = 0.64\lambda_\mathrm{th}`, where :math:`\lambda_\mathrm{th} = \hbar \left/ \sqrt{m k_B T} \right.` is the thermal electron wavelength. If the value after *harm* keyword is -1, then :math:`s_0` mentioned above is used as a harmonic parameter, otherwise the specified value in angstroms is used.

If the *flex_press* keyword is used, then a contribution from the electron widths is added to the total virial and pressure of the system.

This potential is designed to be used with :doc:`atom_style wavepacket <atom_style>` definitions, in order to handle the description of systems with interacting nuclei and explicit electrons.

The following coefficients must be defined for each pair of atoms types via the :doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data file or restart files read by the :doc:`read_data <read_data>` or :doc:`read_restart <read_restart>` commands, or by mixing as described below:

* cutoff (distance units)

For *awpmd/cut*, the cutoff coefficient is optional. If it is not used (as in some of the examples above), the default global value specified in the pair_style command is used.

Currently only Monte Carlo ensemble averaging is supported for awpmd/cut pair style. The Monte-Carlo sampling is given by :doc:`fix mc/wpmd <fix_mc_wpmd>`.


Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for this pair style.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need
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

This pair does not support trajectory integration and works only with :doc:`fix mc/wpmd <fix_mc_wpmd>`.

MPI version has not tested well yet.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

:doc:`pair_style awpmd/cut <pair_awpmd>`

Default
"""""""

These are the defaults for the pair_style keywords: *hartree* for the
initial wave function, *free* for the wave packet width.

----------

.. _Morozov2012:

**(Morozov, 2012)** I.V. Morozov, and I.A. Valuev. Improvement of wave packet molecular dynamics using packet splitting. Contrib. Plasma Phys. 52:140 (2012).

.. _Valuev2015:

**(Valuev, 2015)** I.A. Valuev, and I.V. Morozov. Extension of the wave packet molecular dynamics method towards the accurate quantum simulations of electron dynamics. J. Phys.: Conf. Ser. 653:012153 (2015).

.. _Klakow1994awpmd:

**(Klakow, 1994)** D. Klakow, C. Toepffer, and P.-G. Reinhard, Semiclassical molecular dynamics for strongly coupled coulomb systems, J. Chem. Phys. 101:10766 (1994).

.. _Jakob2007:

**(Jakob, 2007)** B. Jakob, P.-G. Reinhard, C. Toepffer, and G. Zwicknagel. Wave packet simulation of dense hydrogen. Phys. Rev. E. 76:036406 (2007).

.. _Zwicknagel2006:

**(Zwicknagel, 2006)** G. Zwicknagel, and T. Pschiwul. WPMD simulations of a two-component plasma. J. Phys. A. 39:4359 (2006).
