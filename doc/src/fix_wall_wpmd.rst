.. index:: fix wall/wpmd

fix wall/wpmd command
=====================
Syntax
""""""

.. parsed-literal::

   fix ID group-ID wall/wpmd epsilon keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* epsilon = wall potential strength in relative unit :math:`E_0`. (See description)
* zero or more keyword/value pairs may be appended to args

  .. parsed-literal::
         box = length of the box edge L (distance units)
         width_force = take in account width part of force for pressure evaluation
         axes = all or part of *x*, *y* and *z*

Examples
""""""""

.. code-block:: LAMMPS

    fix wall_0 all wall/wpmd 10.0
    fix wall_1 all wall/wpmd 10.0 box 10
    fix wall_2 ion wall/wpmd 10.0 box 10 axes z
    fix wall_3 electron wall/wpmd 10.0 axes x y z

Description
"""""""""""

Restricts the simulation domain with a 3D harmonic potential with center at (0, 0, 0).

.. math::
    E = E_x + E_y + E_z = \epsilon \Big[ \Big( |x| - \frac{L}{2} \Big)^2
        + \Big( |y| - \frac{L}{2} \Big)^2
        + \Big( |z| - \frac{L}{2} \Big)^2 \Big] \\
    E_x = 0,\ \textrm{for}\ -L/2 <= x <= L/2 \\
    E_y = 0,\ \textrm{for}\ -L/2 <= y <= L/2 \\
    E_z = 0,\ \textrm{for}\ -L/2 <= z <= L/2 \\

The potential interacts with particles and wavepackets by generating a force on the particle in
a direction perpendicular to the wall. The potential acts on wavepackets in quantum manner 
affecting both the center and width of the wavepacket. This fix can be used to prevent wavepacket spreading.

The "wall" position is always centered at (0, 0, 0) and its size is determined by *box* keyword setting a length L of the box
edge.

The option *axes* allows to set the axes that are restricted by harmonic potential.

----------


Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

This fix adds the energy of interaction between atoms and all the
specified walls to the global potential energy of the system as part
of :doc:`thermodynamic output <thermo_style>`.

This fix computes a global scalar energy and a global vector with 2 values
, which can be accessed by various :doc:`output commands <Howto_output>`.
The index in vector is:
* 1 energy
* 2 system pressure, evaluated as

.. math::
    P = 1.0 / d \sum_{i=0}^{d} |{F_i}|/L^2,

where i is the coordinate axis index, F_i is the force of the wall, d is the number of enabled walls.

Note that the scalar energy is the sum of interactions with all enabled walls.
If you want the energy on a per-wall basis, you need to use multiple fix wall commands.

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Restrictions
""""""""""""

This fix is part of the AWPMD package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info. The fix oriented to work with atom style :doc:`wavepacket <atom_style>`,
for other atom styles use :doc:`fix wall <fix_wall>`.



Related commands
""""""""""""""""

:doc:`fix wall <fix_wall>`,
:doc:`fix wall/reflect <fix_wall_reflect>`,
:doc:`fix wall/gran <fix_wall_gran>`,
:doc:`fix wall/region <fix_wall_region>`

Default
"""""""

The option defaults box = the minimum length of the simulation region edge,
axes = x y z, wall enabled for all axes.
