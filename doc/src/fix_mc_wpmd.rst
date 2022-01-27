.. index:: fix mc_wpmd

fix mc/wpmd command
================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID mc/wpmd T keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mc/wpmd = style name of this fix command
* T = temperature of the ideal gas reservoir (temperature units)
* zero or more keyword/value pairs may be appended

.. parsed-literal::

     keyword = *ix*, *ex*, *iv*, *ev*, *ew*, *ewp* or *ec*
        *ix* value = vary ions coordinates
        *ex* value = vary wavepacket center coordinates
        *iv* value = vary ions velocity
        *ev* value = vary wavepacket velocity
        *ew* value = vary wavepacket width
        *ewp* value = vary wavepacket width impulse
        *ec* value = vary wavepacket split coefficient

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all mc/wpmd 10000 ix iv ex ev ep epw ec
    fix 1 ions mc/wpmd 5000 ix iv
    fix 1 electrons mc/wpmd 20000 ex ev ep epw

Description
"""""""""""

This fix performs Monte Carlo (MC) moves within the simulation cell or region
for :doc:`wavepacket <atom_style>` atom style. The fix change particle's coordinates,
velocities and wavepacket's width and conjugated impulse at each step. The new configuration
can be accepted or rejected after total energy evaluation. If energy of the new configuration
is grater than before, this configuration can be accepted with probability

.. math::
p = A \exp(-dE/(k_{B}T)).

This fix required arguments that describe operations that will be perform for the particle. On each step
the only one operation will perform.

The *ix* and *ex* options mean that the fix will change particles and packets position on some small values.

The *iv* and *ev* options mean that the fix will change change particles and packets velocities by small values.

The *ew* option is valid only for wavepacket atom style and the fix will change packets width.

The *ewp* option means that impulse related to packet width will be change during the step.

The *ec* option is valid only for wavepacket atom style and split representation of electron. This option means that
fix will change split coefficients for wavepackets.

The amplitude of changes ajusted to provide 50% of accepted steps.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The fix does not writes any information to restarts files.

This fix computes a vector of length 6, which can be accessed by various :doc:`output commands <Howto_output>`.
The vector values are the following global cumulative quantities:

* 1. Accept flag. The step was accepted if this flag is 1 and rejected otherwise.

* 2. The energy of the last accepted configuration.

* 3. The energy of current configuration.

* 4. The number of accepted steps.

* 5. The number of rejected steps.

* 6. The id of the current action over a system. The integer number int range from 0 to k, where k is
a number of options.

Restrictions
""""""""""""

This fix is part of the AWPMD package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info. The fix oriented to work with atom style :doc:`wavepacket <atom_style>`,
the usage with other styles is possible but didn't tested.

Do not set "neigh_modify once yes" or else this fix will never be called. Reneighboring is required.

Can be run in parallel, but some aspects of the MC part will decrease the scale of algorithm.

This fix requires the `thermo  1`.

Use of multiple fix mc/wpmd commands in the same input script can be
problematic due to inconsistency between different fixes.

The neighbor lists are required to re-built every timestep that this fix is
invoked.

Related commands
""""""""""""""""

Default
"""""""

----------
