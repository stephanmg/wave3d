# Wave3d

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e26da086ea2243dbb4399b212885cd4a)](https://app.codacy.com/manual/stephan_5/wave3d?utm_source=github.com&utm_medium=referral&utm_content=stephanmg/wave3d&utm_campaign=Badge_Grade_Dashboard)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)



The repository contains grids (*UGX*) and script (*LUA*) files to initiate a calcium wave in neurons. Additionally job scripts for the local HPC cluster Owlsnest2 at Temple are provided. Required are `SuperLU` and `Parmetis` for parallel computation.

## Prerequisites
Recent versions of [ug4](https://github.com/UG4/ugcore), [neuro collection](https://github.com/NeuroBox3D/neuro_collection), [convection diffusion](https://github.com/UG4/plugin_ConvectionDiffusion), LIMEX, Parmetis and LUA2C are required (The latter three are currently not publicly available but available through the group-internal Quadruped repository).

### Build instructions
First the user needs to enable the plugins above via `cmake -DNEURO_COLLECTION=ON -DLIMEX=ON -DParmetis=ON -DLUA2C=ON -DConvectionDiffusion=ON` and then compile the code via `make`.
An example configuration is displayed [here](config/cmake_setting.txt).

## Running a simulation
Execute in a terminal: 

`ugshell -ex wave3d_revised.lua -grid zhi_y.ugx -dt 1e-8 -endTime 0.0005 -solver GMG -tol 0.02 -ryrDensity 0  -setting ryr`.

For cylinder structure:

`ugshell -ex wave3d_cylinder.lua -grid cylinder_ex.ugx -dt 1e-6 -endTime 0.05 -solver GMG -tol 0.02 -outName .  -ryrDensity 3  -setting ryr -vtk -pstep 0.005`

For y-structure:

`ugshell -ex wave3d_branching.lua -grid y_structure_ex.ugx -dt 1e-6 -endTime 0.05 -solver GMG -tol 0.02 -outName .  -ryrDensity 3  -setting ryr -vtk -pstep 0.005`

For full geometries (`pstep=0` is necessary to write VTK output data in each step):

`ugshell -ex wave3d_revised_sg.lua -grid full_cell_with_soma_and_subsets_assigned.ugx -dt 1e-6 -endTime 0.05 -solver GMG -tol 0.02 -outName . -ryrDensity 3  -setting ryr -vtk -pstep 0.0`

## Notes about LIMEX
-  LIMEX uses C++11 features, so one needs at least the following compiler revisions: [Clang 3.3](https://clang.llvm.org/cxx_status.html) or [GCC 4.8.1](https://gcc.gnu.org/projects/cxx-status.html#cxx11)
- When using LIMEX as the time stepping scheme and concurrently outputting VTK
data one should consider using an older version of the LIMEX plugin since collecting
`pvd` files are not correctly written by: `git checkout 1384a90d0d1f75f582563e71cdf4295f17bfd474`. 

## Example job for Owlsnest
```
#!/bin/bash -l        
#PBS -l walltime=8:00:00,nodes=1:ppn=8,mem=10gb 
#PBS -m abe 
#PBS -M youremail@temple.edu
cd ~/program_directory
module load intel 
module load ompi/intel 
mpirun -np 8 program_name < inputfile > outputfile
```
