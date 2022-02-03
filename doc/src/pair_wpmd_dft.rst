.. index:: pair_style wpmd/dft/cut

pair_style wpmd/dft/cut command
============================

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
Packet Molecular Dynamics (WPMD) method with Hartree approximation[TODO::cite art].

.. math::
    U_{ii}(r) = \frac{C}{r}

    U_{ei}(r, s) = \frac{C}{r} \textrm{Erf}(-\frac{Ar}{s})

    U_{ee}(r, s_1, s_2) = \frac{C}{r} \textrm{Erf}(-\frac{Ar}{f(s_1, s_2)}),

where :math:`A`, :math:`C` --- unit coefficients, :math:`f(s_1, s_2)` is a function[TODO: add full definition].

The pair has only one parameter: `Rc` is the cutoff.

This potential is designed to be used with :doc:`atom_style wavepacket <atom_style>` definitions,
in order to handle the description of systems with interacting nuclei and explicit electrons.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands, or by mixing as described below.
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
This pair is part of the AWPMD package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`
