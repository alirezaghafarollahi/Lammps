# Lammps
 Lammps is a software tool used for molecular dynamics simulations. 

## Steps to perform a lammps simulations 
There are common steps which should be followed in every lammps atomistic simulation. 

- Define units and boundary conditions
- Define geometry
- Define potential parameters including pair/coeff and pair/style
- Define compute, thermo, and dump parameters 
- Define energy minimization scheme, `cg` or `fire`
- Perform the minimization using `minimize etol ftol maxiter maxeval` for 0k and `run N` for finite temperature simulations. `etol` and `ftol` are stopping tolerence for energy and force, respectively. 

Here is an example code for energy minimization at T=0

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

## Defect studies
We can also use Lammps to generate defects such as screw and edge dislocations

Here is an exmple to insert a <111> screw dislocation in a BCC material with periodic boundary conditions

```lammps
clear 
units metal
dimension 3
boundary p s p
atom_style atomic
 
read_data ./structure.inp
include ./potential.inp
 
variable lattice_param equal 3.238
 
change_box all triclinic
 
variable x_spacing equal (sqrt(6.0))
variable y_spacing equal (sqrt(2.0))
variable z_spacing equal (sqrt(3.0)/2)
 
variable x_num_cells equal 5
variable y_num_cells equal 10
variable z_num_cells equal 2
 
lattice bcc ${lattice_param} orient x 1 -1 2 orient y 1  1 0 orient z -1 1 1 origin  0 0 0 spacing ${x_spacing} ${y_spacing} ${z_spacing}
 
variable mid equal -0.1
 
region bot_region block INF INF INF ${mid} INF INF units box
region top_region block INF INF ${mid} INF INF INF units box
 
group bot region bot_region
group top region top_region
############################################################################################
#--- insert dislocation 1 -------------------------------------------------#
############################################################################################
variable Bz equal zlat
variable bz equal ${Bz}
 
variable Xmin_t equal bound(top,xmin)
variable Xmax_t equal bound(top,xmax)

displace_atoms top ramp z -${bz} 0 x ${Xmin_t} ${Xmax_t} units box
 
compute peratom all pe/atom
compute pe      all pe
compute 1 all property/atom x y z
 
variable thermotime equal 100
 
thermo_style custom step temp pe etotal pxx pxy pxz pyy pyz pzz vol fnorm
thermo_modify format float %10.8g
thermo ${thermotime}
 
#Output files
#dump 100 all custom 10000 dump.out.1.* id type x y z c_peratom

##################### minimization #########################
variable iii loop 10
variable tol equal 0.1
 
reset_timestep 0 
label min
 
fix 1 all box/relax x 0.0 z 0.0 xz 0.0
min_style cg
minimize 0.0 1e-4 20000 200000
unfix 1
min_style cg
minimize 0.0 1e-4 20000 200000
 
variable press_x equal abs(pxx)
variable press_z equal abs(pzz)
variable press_xz equal abs(pxz)
variable Fn equal "fnorm"
 
if "${Fn} < 1e-4 && ${press_x} < ${tol} && ${press_z} < ${tol} && ${press_xz} < ${tol}" then "jump SELF exit_label3"
 
next iii
jump SELF min
 
label exit_label3
 
#undump 100
##################### minimization #########################
reset_timestep 0
dump 2 all custom 1 dump.screw.1 id type x y z c_peratom
dump_modify 2 sort id
dump_modify 2 format line "%d %d %.8e %.8e %.8e %.8e"
run 0
undump 2
 
write_data data.screw.1


variable tmp equal "pe"
variable E1 equal ${tmp}
variable str_E1 format E1 "%.10e"
variable tmp delete

print '${str_E1}' append ${energy_file}

quit
```
