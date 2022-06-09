.. index:: fix nve/wpmd

fix nve/wpmd command
=====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nve/wpmd

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/wpmd = style name of this fix command

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/wpmd

Description
"""""""""""

This fix performs the constant NVE integration (constant volume and total energy) to update the positions and velocities of nuclei and electrons in the group for the :doc:`Wave Packet Molecular Dynamics <pair_awpmd>` model. This creates a system trajectory consistent with the microcanonical ensemble.

The operation of this fix is exactly like that described by the :doc:`fix nve <fix_nve>` command, except that the electron radius and radial momentum (width and width momentum of
the wavepackets) are also updated.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the AWPMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix nve <fix_nve>`

Default
"""""""

none
