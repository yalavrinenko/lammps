AWPMD and WPMD-DFT examples
=========================================


This directory contains 8 files with demo input scripts that allow to compute
several properties of atoms, molecules and plasmas via WPMD, AWPMD, and WPMD-DFT methods.

1. `wpmd.H.lammps` -- energy of a single hydrogen atom in Hartree approximation. 
   Here and in 2 - 5 each point in the curve is computed by Monte-Carlo algorithm.
2. `awpmd.H.split.lammps` -- energy of a single hydrogen atom in UHF approximation.
3. `awpmd.H2.lammps` -- potential energy curve for a hydrogen molecule in UHF approximation.
4. `awpmd.H2.split.lammps` -- same as in 3 but with split wavepackets.
5. `wpmc-dft.nvgpu.H2.lammps` -- potential energy curve for hydrogen molecule in WPMD-DFT model,
    required GPU acceleration package. 
6. `awpmd.bulk.lammps` -- energy minimum of bulk hydrogen plasma in AWPMD (UHF approximation) model.
    The configuration with minimal energy is reached in Monte-Carlo algorithm. 
7. `wpmc-dft.nvgpu.bulk.lammps` -- energy minimum of bulk hydrogen plasma in WPMD-DFT model, requires 
    GPU acceleration package. The configuration with minimal energy is reached in Monte-Carlo algorithm.
8. `wpmd-dft.nvgpu.bulk.lammps` -- MD trajectory for equilibrium bulk plasma. There are two
    stages in the simulation:
   1. Equilibration from random configuration via NVT-thermostate.
   2. Time averaging in MD algorithm over equilibrium trajectory
9. `harmonic_trap.lammps` -- single wavepacket in harmonic trap
