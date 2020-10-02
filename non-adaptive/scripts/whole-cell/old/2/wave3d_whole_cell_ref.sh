#!/bin/sh
## PBS setup
#PBS -l walltime=23:59:59
#PBS -V
#PBS -m bea
#PBS -N ONE_REF_WHOLE_CELL
#PBS -q normal
#PBS -l nodes=1:ppn=1

# change to work dir
cd $PBS_O_WORKDIR

# grid and script
SCRIPT=wave3d_whole_cell_ref.lua
GRID=/home/tug41634/new_setups/example_for_markus.ugx

# settings (lowered dt by 1/10, increased num time steps by 10 and vtk write interval by 10 as well)
SOLVER=GS # solver id
DT=1e-7 # delta t
SETTING=all # IP3 and RyR plus Leak
TOL=1e-2
END_TIME=0.100
VTK_FOLDER=vtkWholeCellRef
RYR_DENSITY=1.00
NUM_TIME_STEPS=1000000
VTK_WRITE_INTERVAL=10000
BINARY=ugshell
NP=1
NUM_ER_REFS=0
NUM_PRE_REFS=0
NUM_POST_REFS=0
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
       -writeInterval "$VTK_WRITE_INTERVAL" \
       -numERMRefs "$NUM_ER_REFS" \
       -numPreRefs "$NUM_PRE_REFS" \
       -numPostRefs "$NUM_POST_REFS"
