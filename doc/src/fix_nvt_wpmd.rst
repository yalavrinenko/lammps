.. index:: fix nvt/wpmd

fix nvt/wpmd command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID nvt/wpmd keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *nvt/wpmd*

  .. parsed-literal::

     one or more keyword value pairs may be appended
     keyword = *temp*
       *temp* values = Tstart Tstop Tdamp
         Tstart,Tstop = external temperature at start/end of run
         Tdamp = temperature damping parameter (time units)
       *tchain* value = length of thermostat chain (1 = single thermostat)
       *tloop* value = number of sub-cycles to perform on thermostat
       *nreset* value = reset reference cell every this many timesteps
       *drag* value = drag factor added to barostat/thermostat (0.0 = no drag)

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all nvt/wpmd temp 300.0 300.0 0.1
    fix md all nvt/wpmd temp 10000 10000 $(100*dt)

Description
"""""""""""

These commands perform time integration on Nose-Hoover style
non-Hamiltonian equations of motion for nuclei and electrons in the
group for the :doc:`Wave packet molecular dynamic <pair_awpmd>` model.  The fixes
are designed to generate positions and velocities sampled from the
canonical (nvt) ensembles. The operation of these fixes is exactly like that described by the
:doc:`fix nvt <fix_nh>` command, except that the radius
and radial velocity of electrons are also updated.

.. note::

   Currently, there is no available option for the user to set or
   create temperature distributions that include the radial electronic
   degrees of freedom with the :doc:`velocity <velocity>` command, so the
   the user must allow for these degrees of freedom to equilibrate
   (i.e. equi-partitioning of energy) through time integration.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

See the page for the :doc:`fix nvt, npt, and nph <fix_nh>` commands
for details.

Restrictions
""""""""""""

This fix is part of the AWPMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Other restriction discussed on the page for the :doc:`fix nvt, npt, and nph <fix_nh>` commands also apply.

Related commands
""""""""""""""""

:doc:`fix nvt <fix_nh>`, :doc:`fix nph <fix_nh>`, :doc:`fix npt <fix_nh>`,
:doc:`fix_modify <fix_modify>`, :doc:`run_style <run_style>`

Default
"""""""

See the page for the :doc:`fix nvt, npt, and nph <fix_nh>` commands
for details.

