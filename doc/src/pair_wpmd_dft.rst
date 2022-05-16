.. index:: pair_style wpmd/dft/cut
.. index:: pair_style wpmd/dft/cut/nvgpu

pair_style wpmd/dft/cut command
================================
Accelerator variants: *wpmd/dft/cut/nvgpu*

Syntax
""""""

.. code-block:: LAMMPS

    pair_style wpmd/dft/cut Rc keyword values ...

* Rc = global cutoff, -1 means cutoff of half the shortest box length

.. parsed-literal::

     keyword = *mesh*, *min_cell_size*, *max_distance*, *dynamic*, *force_mesh_bins*
        *mesh* value = *regular* NCells or *adaptive* min_cell_size cell_cutoff
            NCells = number of cells for one direction for regular mesh
            min_cell_size = minimal cell size in fraction of packet width
            cell_cutoff = maximal distance to packet distance for mesh subdivision in packet width
        *dynamic* = on or off
        *force_mesh_cells* value = NCells
            NCells = number of cells for mesh linked with packet
        *gppn* = number of GPU per one node (for accelerator version only)

Examples
""""""""

.. code-block:: LAMMPS

    pair_style wpmd/dft/cut 10.0
    pair_style wpmd/dft/cut 15.0 mesh adaptive 0.7 2.5
    pair_style wpmd/dft/cut 15.0 mesh regular 100
    pair_style wpmd/dft/cut 15.0 mesh adaptive 0.7 2.5 dynamic on force_mesh_cells 11
    pair_style wpmd/dft/cut/nvgpu 15.0 mesh adaptive 0.7 2.5 dynamic on force_mesh_cells 11 gppn 3
    pair_style wpmd/dft/cut/nvgpu 15.0 mesh regular 100 dynamic on force_mesh_cells 11

Description
"""""""""""

This pair style contains an implementation of the :doc:`Wave Packet Molecular Dynamics <pair_wpmd>` method with density functional theory extension (WPMD-DFT), see :ref:`(Lavrinenko, 2021) <wpmddft>`.

The method of WPMD-DFT uses the Hartree approximationfor evaluation of the Coulomb interaction. The total energy is defined as

.. math::
	& E = E_\mathrm{Hartree} + E_\mathrm{a},\\
	& E_\mathrm{a}[n] = \left( T_\mathrm{s}[n] - \sum_{k=1}^{N_\mathrm{e}} T_\mathrm{s}[n_k] \right) + \left( E_\mathrm{xc}[n] - \sum_{k=1}^{N_\mathrm{e}} E_\mathrm{xc}[n_k] \right), \\
	&n_k(r) = \varphi_k(\vec r) \varphi_k^*(\vec r), \qquad
	n(r) = \sum_{k=1}^{N_\mathrm{e}} n_k(r),

where :math:`E_\mathrm{Hartree}` is the Hartree energy (see :doc:`pair_style wpmd/cut <pair_wpmd>`), :math:`E_\mathrm{a}` is an additional exchange-correlation energy given by the kinetic :math:`T_\mathrm{s}` and exchange-correlation :math:`E_\mathrm{xc}` functionals, :math:`n(\vec r)` is the local density, :math:`\varphi_k(\vec{r})` is the wavefunction for :math:`k`-th wavepackets (electron), :math:`N_\mathrm{e}` is the number of electrons. All the functionals are calculated via numerical integration on a 3D space mesh.

In local density approximation with spin (LSDA), the exchange-correlation energy is given by

.. math::
    E_\mathrm{xc}^{\mathrm LSDA}[n_\uparrow,n_\downarrow]= \int\epsilon_\mathrm{xc}(n_\uparrow,n_\downarrow)n (\mathbf{r})\, d\mathbf{r}, \qquad
		n(\mathbf{r}) = n_\uparrow(\mathbf{r}) + n_\downarrow (\mathbf{r}),

The additional kinetic energy of the uniform noninteractive electron gas reads

.. math::
    T_\mathrm{s}[n] = \frac{3}{10}(3\pi^2)^{2/3} \int n(\mathbf{r})^{5/3}\, d\mathbf{r}.

where :math:`n_\uparrow` and :math:`n_\downarrow` are the densities of electrons with the spin up and spin down.

The pair style has the parameters:

*Rc* is the cutoff radius for the Coulomb interaction. For an accurate account of the long-range interaction, it should be greater than the cell size.

The *mesh* keyword defines the type of 3D space mesh. There are two types of mesh implemented: *regular* and *adaptive*. The *regular* option sets the regular mesh with fixed cell size. In this case, the additional parameter *NCells* sets the number of cells for one direction. The total cell size is :math:`\mathrm{NCells}`. The *adaptive* option sets the adaptive mesh with variable cell size that depends on the local gradient of the electron density. In this case, the cell width will be greater or equal *min_cell_size*. The parameter *cell_cutoff* defines the maximum distance between the cell center and the wavepacket. The adaptive mesh refinement algorithm improves the simulation performance as the total number of cells decreases.

The *dynamic* keyword enables the force calculation from exchange-correlation interaction. The forces are calculated by the numerical integration at the wavepacket position

.. math::

    \frac{\partial E_{\mathrm{a}}}{\partial q} =
    \left(\frac{\partial T_\mathrm{s}[n]}{\partial n} + \frac{\partial E_\mathrm{xc}[n]}{\partial n}\right)
    \frac{\partial n}{\partial q},

where :math:`q` is one of the wavepacket parameters (position, momentum, width, and momentum of the width).

The *force_mesh_cell* keyword sets the number of cells used for the calculation of force acting on a single wavepacket in one direction. The total number of cells is a cube of *force_mesh_cell*.

The *gppn* keyword sets the number of GPUs per MPI process. This is required for the correct task scattering to multiple GPUs.

This potential inherits all properties of :doc:`pair_style wpmd/cut <pair_wpmd>`.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options are not relevant for this pair style. For minimization of the system's energy with `wpmd/cut` pair should use a :doc:`fix mc/wpmd <fix_mc_wpmd>`.

This pair style writes its information to :doc:`binary restart files <restart>`, so pair_style and pair_coeff commands do not need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command. It does not support the *inner*, *middle*, and *outer* keywords.

----------

Restrictions
""""""""""""
This pair is part of the AWPMD package.  It is only enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

This pair is work only with *real* and *electron* units due to energy conversion units.

This pair required :doc:`wavepacket <atom_style>` or :doc:`electron <atom_style>` atom style.

The time evolution can be calculated by one of :doc:`fix nvt/wpmd <fix_nh_wpmd>`, :doc:`fix nph/wpmd <fix_nh_wpmd>`, :doc:`fix npt/wpmd <fix_nh_wpmd>` or :doc:`fix nve/wpmd <fix_nve_wpmd>`. The Monte-Carlo sampling can be performed with :doc:`fix mc/wpmd <fix_mc_wpmd>`.

GPU acceleration is available for NVidia GPUs only.


Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

:doc:`pair_style wpmd/cut <pair_wpmd>`

.. _wpmddft:


Default
"""""""

By default the *mesh* is *regular* with the 50 cells in each direction. The *dynamic* is *off*.

----------

**(Lavrinenko, 2021)** Lavrinenko, Yaroslav, et al. "Equilibrium properties of warm dense deuterium
calculated by the wave packet molecular dynamics and density functional theory method."
Physical Review E 104.4 (2021): 045304.
