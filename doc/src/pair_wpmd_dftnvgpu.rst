.. index:: pair_style wpmd/dft-nvgpu/cut

pair_style wpmd/dft-nvgpu/cut command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style wpmd/dft-nvgpu/cut Rc keyword values ...

* Rc = global cutoff, -1 means cutoff of half the shortest box length

.. parsed-literal::

     keyword = *mesh*, *min_cell_size*, *max_distance*, *dynamic*, *force_mesh_bins*, *gppn*
        *mesh* value = *regular* NCells or *adaptive* min_cell_size cell_cutoff
            NCells = number of cells for one direction for regular mesh
            min_cell_size = minimal cell size in fraction of packet width
            cell_cutoff = maximal distance to packet distance for mesh subdivision in packet width
        *dynamic* = on or off
        *force_mesh_cells* value = NCells
            NCells = number of cells for mesh linked with packet
        *gppn* = number of GPU per one node
Examples
""""""""

.. code-block:: LAMMPS

    pair_style wpmd/dft/cut 10.0
    pair_style wpmd/dft/cut 15.0 mesh adaptive 0.7 2.5
    pair_style wpmd/dft/cut 15.0 mesh regular 100
    pair_style wpmd/dft/cut 15.0 mesh adaptive 0.7 2.5 dynamic on force_mesh_cells 11

Description
"""""""""""

This pair style contains an implementation of the Wave
Packet Molecular Dynamics (WPMD) method with density functional theory extension :ref:`(wpmddft) <wpmddft>`.

WPMD-DFT uses Hartree approximation (same as in :doc:`pair wpmd/cut <wpmd/cut>`) for evaluation of Coulomb interaction. The additional exchange-correlation
energy evaluates via numerical integration of xc-functionals on 3d mesh.

The numerical integration perform on GPU. At this moment there are only nvidia gpus are supported.

.. math::
    E = E_\mathrm{Hartree} + E_\mathrm{a}

    E_\mathrm{a}[n] = (T_\mathrm{s}[n] - \sum_{i} T_\mathrm{s}[n_i]) +
                    (E_\mathrm{XC}[n] - \sum_{i} E_\mathrm{XC}[n_i])

    n = n(r) = \sum_{k=1}^{N_\mathrm{e}} \varphi(\vec{r}) \varphi^*(\vec{r})

    n_i = n_i(r) = \varphi(\vec r) \varphi^*(\vec r)

where :math:`N_\mathrm{e}` --- number of wavepackets.

The exchange-correlation energy evaluated in local density approximation with spin as:

.. math::
    E_{\mathrm XC}^{\mathrm LSDA}[n_\uparrow,n_\downarrow]=
        \int\epsilon_{\mathrm XC}(n_\uparrow,n_\downarrow)n (\mathbf{r})\, d\mathbf{r},

        n(\mathbf{r}) = n_\uparrow(\mathbf{r}) + n_\downarrow (\mathbf{r}),

Additional kinetic energy of uniform noninteractive electon gas is:

.. math::
    T_\mathrm{s}[n] = \frac{3}{10}(3\pi^2)^{2/3} \int n(\mathbf{r})^{5/3}\, d\mathbf{r}.

The pair has several parameters similar to :doc:`pair_style <wpmd/dft/cut>`:

* The *Rc* is the cutoff radius for Coulomb interaction. Due to accurate account of long range interaction
should be grater then cell size.

* The *mesh* keyword is set up type of 3d space mesh for numerical integration. There are two types of meshes
ware implemented: *regular* and *adaptive*. The *regular* option sets the regular mesh with fixed cell size.
Additional parameter *NCells* set the number of cells for one direction. The total cell size is
:math:`\mathrm{NCells}`. The *adaptive* option sets the adaptive mesh with variable cell size that depend
on the gradient of electron density. The cell width will be grater or equal *min_cell_size*. The parameter
*cell_cutoff* define the maximum distance from cell center to packet center. The adaptive mesh refinement
algorithm increase a performance of simulation due to decreasing of cells number.

* The *dynamic* keyword is enable a force calculation from exchange-correlation interaction. Forces calculates
by numerical integration over a mesh linked to wavepacket.
.. math::
      \frac{\partial E_{\mathrm{a}}}{\partial q} =
    \left(\frac{\partial T_\mathrm{s}[n]}{\partial n} + \frac{\partial E_\mathrm{XC}[n]}{\partial n}\right)
    \frac{\partial n}{\partial q}.

* The *force_mesh_cell* keyword is set the number of cells for force calculation per packet in one direction.
The total number of cells is a cube of *force_mesh_cell*.

* The *gppn* keyword is set the number of gpu per one node. This is required for correct task scattering to gpu.

This potential inherit all properties of :doc:`pair wpmd/cut <wpmd/cut>`.

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The :doc:`pair_modify <pair_modify>` mix, shift, table, and tail options
are not relevant for this pair style. For minimization of the system's energy with `wpmd/cut` pair
should use a :doc:`fix mc/wpmd <fix mc_wpmd>`.

This pair style writes its information to :doc:`binary restart files <restart>`,
so pair_style and pair_coeff commands do not need
to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.
----------

Restrictions
""""""""""""
This pair is part of the WPMD-NVGPU-DFT package.  It is only enabled if LAMMPS was
built with that package. This pair requires a nvidia cuda version 9.0 or higher and c++14.
See the :doc:`Build package <Build_package>` doc page for more info.

This pair is work only with *real* unit due to energy conversion units.
Default
"""""""
By default the *mesh* is *regular* with size per one axe is equal to 50. The *dynamic* is *off*.
The value of *gppn* is 1.

----------

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`
:doc:`pair_style <wpmd/dft/cut>`

.. _wpmddft:

**(wpmddft)** Lavrinenko, Yaroslav, et al. "Equilibrium properties of warm dense deuterium
calculated by the wave packet molecular dynamics and density functional theory method."
Physical Review E 104.4 (2021): 045304.
