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
        *ec_re* value = vary wavepacket split real coefficients (for multi-wavepacket electrons)
        *ec_im* value = vary wavepacket split imag coefficients (for multi-wavepacket electrons)

Examples
""""""""

.. code-block:: LAMMPS

    fix 1 all mc/wpmd 10000 ix iv ex ev ep epw ec
    fix 1 ions mc/wpmd 5000 ix iv
    fix 1 electrons mc/wpmd 20000 ex ev ep epw

Description
"""""""""""

This fix performs the Monte-Carlo (MC) moves within the simulation cell or region for the :doc:`wavepacket <atom_style>` atom style (see :doc:`pair_style wpmd/cut <pair_wpmd>`). Electrons in this atom style may be represented by single or multiple Gaussian wavepackets. For models with one wavepacket per electron, the wavepacket center and width are identical with electron radius and radius momentum correspondingly. The fix changes particles' coordinates, velocities, and wavepackets' width and conjugated momenta at each step forming a trial particle configuration. The new configuration can be accepted or rejected after the evaluation of the total energy. If the energy of the new configuration is greater than the energy of the previously accepted configuration, the current configuration can be accepted with the probability

.. math::
    p = A \exp(-dE/(k_{B}T)).

This fix requires arguments that describe operations on particles while forming a trial configuration. At each step only one of the following operations is performed:

For *ix* and *ex* options, the fix varies the wavepacket (electron) and other particle positions.

For *iv* and *ev* options, the fix varies the wavepacket and other particle velocities.

For *ew* option, the fix varies the wavepacket widths (applicable for wavepacket atom style only).

For *ewp* option, the fix varies the momenta related to the wavepacket widths (applicable for wavepacket atom style only).

For *ec* option, the fix varies the normalizing coefficients of the wavepackets within an electron (applicable for multiple wavepackets per electron only).

The amplitude of adjusted is continuously ajusted to provide 50% of accepted steps in average.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The fix does not write any information to the restart files.

This fix computes a vector of length 6, which can be accessed by various :doc:`output commands <Howto_output>`. The vector values are the following global cumulative quantities:

    1. Accept flag. The step was accepted if this flag is 1 and rejected otherwise.
    2. The energy of the last accepted configuration.
    3. The energy of the current configuration.
    4. The number of accepted steps.
    5. The number of rejected steps.
    6. The id of the current operation performed over a system in the current step. The integer number int ranges from 0 to :math:`k`, where :math:`k` is the number of possible operations.

Restrictions
""""""""""""

This fix is a part of the AWPMD package. It is only enabled if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info. The fix is designed to work with atom styles :doc:`wavepacket or electron <atom_style>`. The usage with other styles is possible but has not been tested.

Do not set "neigh_modify once yes" or else this fix will never be called. Reneighboring is required.

Can be run in parallel, but some aspects of the MC part will decrease the parallel scaling of algorithm.

This fix requires the `thermo  1` command.

Use of multiple fix mc/wpmd commands in the same input script can be problematic due to inconsistency between different fixes.

The neighbor lists are required to re-built every timestep that this fix is invoked.

Related commands
""""""""""""""""

:doc:`fix nvt/wpmd <fix_nh_wpmd>`, :doc:`fix npt/wpmd <fix_nh_wpmd>`, :doc:`fix nph/wpmd <fix_nh_wpmd>`

Default
"""""""

none