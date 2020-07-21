# Wave3d

The repository contains grids (*UGX*) and script (*LUA*) files to initiate a calcium wave in neurons.

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

## Notes about LIMEX
- LIMEX uses C++11 features, so one needs at least the following compiler revisions: 
[Clang 3.3](https://clang.llvm.org/cxx_status.html) or [GCC 4.8.1](https://gcc.gnu.org/projects/cxx-status.html#cxx11)
- When using LIMEX as the time stepping scheme and concurrently outputting VTK
data one needs to use an older version of the LIMEX plugin before running the
simulation by: `git checkout 1384a90d0d1f75f582563e71cdf4295f17bfd474`. This
is a bug and should be addressed soon by the Limex plugin devlopers hopefully.
