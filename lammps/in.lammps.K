units lj
atom_style full

lattice fcc 0.8
region box block 0 6 0 6 0 6
create_box 2 box
create_atoms 1 box
mass 1 1.0
set group all type/subset 2 64 2324349
mass 2 1.0

variable tstep equal 0.001
timestep ${tstep}
variable Tdamp equal 100*${tstep}

variable tinit equal 10.0
variable tfinal equal 5.0

velocity all create ${tinit} 87287

group A type 1
group B type 2
neighbor 0.3 bin
neigh_modify every 5 delay 0 check yes

pair_style lj/cut 5
pair_coeff 1 1 1 1 5
pair_coeff 1 2 1.4142 1.5 5
pair_coeff 2 2 2 2 5
pair_modify shift yes

thermo_style custom step temp pe ke evdwl epair etotal press
thermo_modify flush yes
thermo 100
dump trj all custom 20001 init.lammpstrj id type mol x y z
dump_modify trj format line "%d %d %d %.15f %.15f %.15f"

fix integra all nve/limit 0.05
run 20000
unfix integra
undump trj

fix rescala all temp/rescale 200 ${tinit} ${tinit} 0.02 0.5
fix integra all nve
run 10000
unfix rescala
run 100000
unfix integra

thermo 1

compute myRDF all rdf 500 1 1 1 2 2 2
compute r2t all msd

fix msd1 all ave/time 1 2000 2000 c_r2t[4] file msd_T_${tfinal}.profile 
fix crdf all ave/time 1 10000 10000 c_myRDF[*] file grT_${tfinal}.rdf mode vector

fix integra all nvt temp ${tfinal} ${tfinal} ${Tdamp}
dump trj all custom 10 trj_${tfinal}.lammpstrj id type mol x y z
dump_modify trj format line "%d %d %d %.15f %.15f %.15f"
### NVT ###
run 100000
