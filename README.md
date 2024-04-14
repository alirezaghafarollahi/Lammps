# Lammps
 Lammps is a software tool used for molecular dynamics simulations. 

## Steps to perform simulations 
There are common steps which should be followd in every atomistic simulation. 

- Define units and boundary conditions
- Define geometry depending on the crystal structure
- Define potential parameters including pair/coeff and pair/style
- Define compute, thermo, and dump styles 
- Define energy minimization scheme, "cg" or "fire", and parameters such as energy and force criteria. Good values for energy criteri is 0 and for force is 1e-5.
- Perform the minimization using minimize for 0k and run for finite temperature simulations.

Here is an example code 

``` lammps
clear
units metal
dimension 3
boundary p p p
atom_style atomic
variable a equal 4.0
lattice fcc ${a}

region 1 block -3.0 3.0 -3.0 3.0 -4.0 4.0 untis lattice
create_box 1 box
create_atoms 1 box
pair_style eam/alloy
pair_coeff * * Al.eam.alloy Al # Al.eam.alloy is the name of the potential file

compute peratom all pe/atom
compute pe all pe

thermo_style custom step temp pe vol press pxx pyy pzz fnorm
themo_modify format float %5.5g
variable thermo_time equal 100 # print thermo every N times
thermo ${thermo_time}

variable dump_id equal 1
variable N equal 100 #dump on timesteps which are multiples of N
variable dump_out_name str "dump.out"
dump ${dump_id} all custom ${N} ${dump_out_name}.* id type x y z c_peratom

fix 1 all box/relax iso 0
min_style cg
minimize 0 1e-5 20000 20000
unfix 1
```
## Material properties
Lammps can be used to compute the material properties of metals such as:
 - Lattice constant
 - Cohesive energy
 - Surface energy

### Lattice constant
Lattice constant is commonly computed at T=0k. 
-To compute lattice constant we must use "fix 1 all box/relax iso 0" to let the pressure relax in all directions, where we have periodic boundary conditions in all directions. 
### Cohesive energy
-Same as lattice constant, cohesive energy is computed at T=0 K.
-Cohesive energy is the total potential divided by the toral number of atoms.
### Surface energy
- Surface energy is computed by the subtraction of two potnetial energies. One with periodic boundary condition and the other with free coundary conditions divided by the area.

Here is an example code snippet to compute surface energy: 

``` lammps
clear
units metal
dimension 3
boundary p p p
atom_style atomic
variable a equal 4.0
lattice fcc ${a}

region 1 block -3.0 3.0 -3.0 3.0 -4.0 4.0 untis lattice
create_box 1 box
create_atoms 1 box
pair_style eam/alloy
pair_coeff * * Al.eam.alloy Al # Al.eam.alloy is the name of the potential file

compute peratom all pe/atom
compute pe all pe

thermo_style custom step temp pe vol press pxx pyy pzz fnorm
themo_modify format float %5.5g
variable thermo_time equal 100 # print thermo every N times
thermo ${thermo_time}

variable dump_id equal 1
variable N equal 100 #dump on timesteps which are multiples of N
variable dump_out_name str "dump.out"
dump ${dump_id} all custom ${N} ${dump_out_name}.* id type x y z c_peratom

fix 1 all box/relax iso 0
min_style cg
minimize 0 1e-5 20000 20000
unfix 1

variable tmp equal pe
variable E1 equal ${tmp}
variable tmp delete

change_box all boundary p s p # If we want to compute the surface energy of y-plane
min_style cg
minimize 0 1e-5 20000 20000

variable tmp equal pe
variable E2 equal ${tmp}
variable tmp delete

variable Area equal "lx"*"lz"

variable surface_energy equal (${E2}-${E1})/(2*${Area}) $ 2 is because we have two surfaces.

print "surface energy is: ${surface_energy} "
``` 
