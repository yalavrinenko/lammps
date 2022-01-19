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

This pair style contains an implementation of the Wave
Packet Molecular Dynamics (WPMD) method with Hartree approximation[TODO::cite art].

.. math::
    U_{ii} = \frac{1}{{4\pi \varepsilon _0 }} \sum_{i<j} \frac{Z_i Z_j}{R_{ij}} \\
    U_{ei} = -\frac{1}{{4\pi \varepsilon _0 }} \sum_{i,j} \frac{Z_j e^2}{R_{ij}} \textrm{erf} \Big( -\frac{\sqrt{3} R_{ij}}{\sqrt{2} s_i} \Big) \\
    U_{ee} = \frac{1}{{4\pi \varepsilon _0 }} \sum_{i,j} \frac{e^2}{R_{ij}} \textrm{erf} \Big( -\frac{\sqrt{3} r_{ij}}{\sqrt{2(s_i^2 + s_j^2)}} \Big),

where :math:`Z` --- is an ion charge, :math:`R_{ij}` --- distance between ions or ion and electron,
:math:`r_{ij}` --- distance between electrons.

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

The system evolution in time should be perform by fix with :doc:`*/wpmd <fix_nh_wpmd>` suffix.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

:doc:`pair_style awpmd/cut <pair_awpmd>`

:doc:`pair_style eff/cut <pair_eff>`