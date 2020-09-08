#!/bin/sh
#PBS -l walltime=06:00:00
#PBS -V
#PBS -m bea
#PBS -N wave3d_whole_Cell
#PBS -q normal
#PBS -l nodes=1:ppn=1
cd $PBS_O_WORKDIR

# grid and script
SCRIPT=markus.lua
GRID=/home/tug41634/ExampleSubsets.ugx

# settings
SOLVER=GS # solver id
DT=1e-5 # delta t
SETTING=all # IP3 and Ryr plus Leak
TOL=1e-2
END_TIME=0.100
VTK_FOLDER=vtkWholeCell
RYR_DENSITY=1.00
NUM_TIME_STEPS=10000
VTK_WRITE_INTERVAL=100
mpirun -np 1 ugshell \
       -ex wave3d_whole_cell.lua \
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
