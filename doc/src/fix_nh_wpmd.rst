.. index:: fix nvt/wpmd
.. index:: fix npt/wpmd
.. index:: fix nph/wpmd

fix nvt/wpmd command
===================

fix npt/wpmd command
===================

fix nph/wpmd command
===================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID style_name keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style_name = *nvt/wpmd* or *npt/wpmd* or *nph/wpmd*

  .. parsed-literal::

     one or more keyword value pairs may be appended
     keyword = *temp* or *iso* or *aniso* or *tri* or *x* or *y* or *z* or *xy* or *yz* or *xz* or *couple* or *tchain* or *pchain* or *mtk* or *tloop* or *ploop* or *nreset* or *drag* or *dilate*
       *temp* values = Tstart Tstop Tdamp
         Tstart,Tstop = external temperature at start/end of run
         Tdamp = temperature damping parameter (time units)
       *iso* or *aniso* or *tri* values = Pstart Pstop Pdamp
         Pstart,Pstop = scalar external pressure at start/end of run (pressure units)
         Pdamp = pressure damping parameter (time units)
       *x* or *y* or *z* or *xy* or *yz* or *xz* values = Pstart Pstop Pdamp
         Pstart,Pstop = external stress tensor component at start/end of run (pressure units)
         Pdamp = stress damping parameter (time units)
       *couple* = *none* or *xyz* or *xy* or *yz* or *xz*
       *tchain* value = length of thermostat chain (1 = single thermostat)
       *pchain* values = length of thermostat chain on barostat (0 = no thermostat)
       *mtk* value = *yes* or *no* = add in MTK adjustment term or not
       *tloop* value = number of sub-cycles to perform on thermostat
       *ploop* value = number of sub-cycles to perform on barostat thermostat
       *nreset* value = reset reference cell every this many timesteps
       *drag* value = drag factor added to barostat/thermostat (0.0 = no drag)
       *dilate* value = *all* or *partial*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nvt/wpmd temp 300.0 300.0 0.1
   fix 1 part npt/wpmd temp 300.0 300.0 0.1 iso 0.0 0.0 1.0
   fix 2 part npt/wpmd temp 300.0 300.0 0.1 tri 5.0 5.0 1.0
   fix 2 ice nph/wpmd x 1.0 1.0 0.5 y 2.0 2.0 0.5 z 3.0 3.0 0.5 yz 0.1 0.1 0.5 xz 0.2 0.2 0.5 xy 0.3 0.3 0.5 nreset 1000

Description
"""""""""""

These commands perform time integration on Nose-Hoover style
non-Hamiltonian equations of motion for nuclei and electrons in the
group for the :doc:`wavepacket molecular dynamic <pair_wpmd>` model.  The fixes
are designed to generate positions and velocities sampled from the
canonical (nvt), isothermal-isobaric (npt), and isenthalpic (nph)
ensembles.  This is achieved by adding some dynamic variables which
are coupled to the particle velocities (thermostatting) and simulation
domain dimensions (barostatting).  In addition to basic thermostatting
and barostatting, these fixes can also create a chain of thermostats
coupled to the particle thermostat, and another chain of thermostats
coupled to the barostat variables. The barostat can be coupled to the
overall box volume, or to individual dimensions, including the *xy*,
*xz* and *yz* tilt dimensions. The external pressure of the barostat
can be specified as either a scalar pressure (isobaric ensemble) or as
components of a symmetric stress tensor (constant stress ensemble).
When used correctly, the time-averaged temperature and stress tensor
of the particles will match the target values specified by
Tstart/Tstop and Pstart/Pstop.

The operation of these fixes is exactly like that described by the
:doc:`fix nvt, npt, and nph <fix_nh>` commands, except that the radius
and radial velocity of electrons are also updated.  Likewise the
temperature and pressure calculated by the fix, using the computes it
creates (as discussed in the :doc:`fix nvt, npt, and nph <fix_nh>`
doc page), are performed with computes that include the WPMD contribution
to the temperature or kinetic energy from the electron radial velocity.

.. note::

   currently, there is no available option for the user to set or
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

The keyword defaults are tchain = 3, pchain = 3, mtk = yes, tloop =
ploop = 1, nreset = 0, drag = 0.0, dilate = all, and couple = none.
