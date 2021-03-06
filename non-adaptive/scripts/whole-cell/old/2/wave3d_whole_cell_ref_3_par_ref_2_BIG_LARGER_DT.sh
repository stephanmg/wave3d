#!/bin/sh
## PBS setup
#PBS -l walltime=47:59:59
#PBS -V
#PBS -m bea
#PBS -N BIG_THREE_REFS_WHOLE_CELL_PAR_REF_TWO_BIG_LARGER_DT
#PBS -q normal
#PBS -l nodes=12:ppn=28

# change to work dir
cd $PBS_O_WORKDIR

# grid and script
SCRIPT=wave3d_whole_cell_ref.lua
GRID=/home/tug41634/new_setups/example_markus_3rd_refinement.ugx

# settings
SOLVER=GS # solver id
DT=1e-5 # delta t
SETTING=all # IP3 and Ryr plus Leak
TOL=1e-2
END_TIME=0.100
VTK_FOLDER=vtkWholeCellRefThreeParRefTwoBigLargerDt
RYR_DENSITY=1.00
NUM_TIME_STEPS=10000000
VTK_WRITE_INTERVAL=100000
BINARY=ugshell
NP=336
NUM_ER_REFS=0
NUM_PRE_REFS=0
NUM_POST_REFS=2
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
