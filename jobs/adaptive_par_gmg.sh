#!/bin/sh
## PBS setup
#PBS -l walltime=23:59:59
#PBS -V
#PBS -m bea
#PBS -N NEW_ADAPTIVE_GMG
#PBS -q normal
#PBS -l nodes=1:ppn=28

# change to work dir
cd $PBS_O_WORKDIR

# driver script
SCRIPT=wave3d_whole_cell_adaptive.lua
# grids
GRID=/home/tug41634/new_adaptive_setup/synapse_simulation_benchmark_3d_input.ugx
# output
VTK_FOLDER=adaptive_gmg_par
# simulation settings and aprameters
SOLVER=GMG # solver id
DT=1e-6 # delta t 
SETTING=all # IP3 and Ryr plus Leak
END_TIME=1.00  # end time
RYR_DENSITY=1.00 # ryanodine receptors
TOL=1e-12

# execution parameters
BINARY=ugshell
NP=28 # 28 worked
# execute simulation
mpirun -np "$NP" "$BINARY" \
       -ex "$SCRIPT"  \
       -grid "$GRID"  \
       -endTime "$END_TIME" \
       -vtkFolder "$VTK_FOLDER" \
       -dt "$DT" \
       -setting "$SETTING" \
       -solver "$SOLVER"  \
       -tol "$TOL" \
       -pstep 0.0 \
       -writeInterval 1 \
       -numERMRefs 0 \
       -numPreRefs 0 \
       -numTimeSteps 10000
