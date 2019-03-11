units           lj
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 2 0 2 0 2
create_box      1 box
create_atoms    1 random 128 123123 box
mass            1 1.0

velocity        all create 3.0 87287

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5

compute         energies all pair lj/cut
compute         kinetik all ke

variable        mc_eng equal c_energies+c_kinetik
fix             mc all wpmc/awpmd mc_eng 1

variable                accept equal f_mc[1]
variable                mc_energy equal f_mc[2]
variable                step_energy equal f_mc[3]

fix 			ave all ave/time 1 1 1 c_energies c_kinetik v_mc_eng v_accept c_energies[*] file ave.dat off 1

thermo          1
thermo_style    custom step etotal pe ke temp v_accept v_mc_energy v_step_energy
dump h5md1 		all h5md 100 dump_h5md.h5 position velocity

run 10000000
