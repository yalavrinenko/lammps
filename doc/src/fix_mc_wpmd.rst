.. index:: fix mc/wpmd

fix mc/wpmd command
====================

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
        *ix* value = vary ion coordinates
        *ex* value = vary wavepacket center coordinates (electron positions)
        *iv* value = vary ion velocities
        *ev* value = vary wavepacket velocities
        *ew* value = vary wavepacket widths (electron radii)
        *ewp* value = vary wavepacket width momenta  (electron radius momenta)
        *ec* value = vary wavepacket split coefficients (for multi-wavepacket electrons)

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all mc/wpmd 10000 ix iv ex ev ep epw ec
    fix 1 ions mc/wpmd 5000 ix iv
    fix 1 electrons mc/wpmd 20000 ex ev ep epw

Description
"""""""""""

This fix performs Monte Carlo (MC) moves within the simulation cell or region
for :doc:`wavepacket <atom_style>` atom style. Electrons in this atom style may be represented by a sinle or multiple wavepackets. 
For models with one wavepacket per electron the wavepacket center and width are identical with elcetron radius and radius momentum correspondingly.
The fix changes particles' coordinates,
velocities and wavepackets' width and conjugated momenta at each step forming a trial particle configuration. The new configuration
can be accepted or rejected after total energy evaluation. If the energy of the new configuration
is grater than the enery of the previously accepted configuration, the current configuration can be accepted with probability

.. math::
    p = A \exp(-dE/(k_{B}T)).

This fix requires arguments that describe operations that will be performed for the particle while forming a trial confiuration. At each step
only one of the following operations are performed.

The *ix* and *ex* options mean that the fix will shift particles and electron (wavepacket) positions by some small values.

The *iv* and *ev* options mean that the fix will change particles and packets velocities by small values.

The *ew* option is valid only for wavepacket atom style and the fix will change packets widths.

The *ewp* option means that momenta related to the packet widths will be changed during the step.

The *ec* option is valid only for wavepacket atom style and the split representation of electrons. This option means that
fix will change split coefficients for wavepackets.

The amplitude of changes is continuously ajusted to provide 50% of accepted steps in average.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The fix does not write any information to restarts files.

This fix computes a vector of length 6, which can be accessed by various :doc:`output commands <Howto_output>`.
The vector values are the following global cumulative quantities:

    * 1. Accept flag. The step was accepted if this flag is 1 and rejected otherwise.
    * 2. The energy of the last accepted configuration.
    * 3. The energy of the current configuration.
    * 4. The number of accepted steps.
    * 5. The number of rejected steps.
    * 6. The id of the current operation performed over a system in the current step. The integer number int ranges from 0 to k, where k is the number of possible operations.

Restrictions
""""""""""""

This fix is part of the AWPMD package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
doc page for more info. The fix is designed to work with atom styles :doc:`wavepacket or electron <atom_style>`.
The usage with other styles is possible but has not been tested.

Do not set "neigh_modify once yes" or else this fix will never be called. Reneighboring is required.

Can be run in parallel, but some aspects of the MC part will decrease the parallel scaling of algorithm.

This fix requires the `thermo  1` command.

Use of multiple fix mc/wpmd commands in the same input script can be
problematic due to inconsistency between different fixes.

The neighbor lists are required to re-built every timestep that this fix is
invoked.

Related commands
""""""""""""""""

Default
"""""""

