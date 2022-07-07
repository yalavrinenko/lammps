AWPMD and WPMD-DFT examples
=========================================

This directory contains input scripts that allow to compute properties of atoms, molecules and plasmas using WPMD, AWPMD, and WPMD-DFT methods.

`in.wpmd.H.mc` calculates the ground state energy of the hydrogen atom using the original Wave Packet Monte Carlo algorithm with a single Gaussian per electron. 

`in.awpmd.H.split.mc` calculates the ground state energy of the hydrogen atom using the Wave Packet Monte Carlo algorithm with multiple Gaussians per electron.

`in.awpmd.H2.mc` calculates the ground state energy of H2 molecule depending on the interatomic distance using the Antisymmetrized Wave Packet Monte Carlo (UHF approximation) with a single Gaussian per electron.

`in.awpmd.H2.split.mc` calculates the ground state energy of H2 molecule depending on the interatomic distance using the Antisymmetrized Wave Packet Monte Carlo (UHF approximation) with multiple Gaussians per electron.

`in.awpmd.bulk.mc` calculates the energy minimum of a bulk hydrogen plasma using the Antisymmetrized Wave Packet Monte Carlo (UHF approximation) with a single Gaussian per electron.

`in.awpmd.bulk.split.mc` calculates the energy minimum of a bulk hydrogen plasma using the Antisymmetrized Wave Packet Monte Carlo (UHF approximation) with multiple Gaussians per electron.

`in.wpmd-dft.H2.mc` calculates the ground state energy of H2 molecule depending on the interatomic distance using the WPMC-DFT algorithm.

`in.wpmd-dft.nvgpu.H2.mc` calculates the ground state energy of H2 molecule depending on the interatomic distance using the WPMC-DFT algorithm. GPU acceleration package is required.

`in.wpmd-dft.bulk.mc` calculates the energy minimum of a bulk hydrogen plasma using the WPMC-DFT algorithm.
    
`in.wpmd-dft.bulk.md` calculates the MD trajectory for equilibrium bulk plasma. There simulation consists of two stages:
   1. Equilibration from a random configuration using NVT thermostat.
   2. Calculation of the equilibrium trajectory with time averaging of the thermodynamic properties

`in.wpmd-dft.nvgpu.bulk.mc` calculates the energy minimum of a bulk hydrogen plasma using the WPMC-DFT algorithm. GPU acceleration package is required.

`in.wpmd-dft.nvgpu.bulk.md` calculates the MD trajectory for equilibrium bulk plasma. GPU acceleration package is required. 

`in.wpmd.harmonic_trap.mc` single wavepacket in a harmonic trap
