#!/bin/sh
## PBS setup
#PBS -l walltime=01:59:59
#PBS -V
#PBS -m bea
#PBS -N NO_REF_PARALLEL
#PBS -q normal
#PBS -l nodes=1:ppn=8

# change to work dir
cd $PBS_O_WORKDIR

# grid and script
SCRIPT=wave3d_whole_cell_par.lua
GRID=/home/tug41634/new_setups/example_for_markus_not_refined.ugx

# settings
SOLVER=GS # solver id
DT=1e-6 # delta t
SETTING=all # IP3 and Ryr plus Leak
TOL=1e-2
END_TIME=0.100
VTK_FOLDER=vtkWholeCellPar
RYR_DENSITY=1.00
NUM_TIME_STEPS=100000
VTK_WRITE_INTERVAL=1000
BINARY=ugshell
NP=8
mpirun -np "$NP" "$BINARY" \
       -ex "$SCRIPT"  \
       -grid "$GRID"  \
       -endTime "$END_TIME" \
       -tol "$TOL" \
       -vtk  \
       -vtkFolder "$VTK_FOLDER" \
       -dt "$DT" \
       -setting "$SETTING" \
       -ryrDensity "$RYR_DENSITY" \
       -solver "$SOLVER"  \
       -numTimeSteps "$NUM_TIME_STEPS"  \
       -writeInterval "$VTK_WRITE_INTERVAL"
