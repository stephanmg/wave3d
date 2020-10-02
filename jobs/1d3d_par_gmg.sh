#!/bin/sh
## PBS setup
#PBS -l walltime=05:59:59
#PBS -V
#PBS -m bea
#PBS -N PAR_GMG_2_REGIONS
#PBS -q normal
#PBS -l nodes=2:ppn=28

# change to work dir
cd $PBS_O_WORKDIR

# driver script
SCRIPT=1d3d_singleCell_synDistr_new_ball.lua
# grids
GRID_3D=/home/tug41634/two_regions/synapse_simulation_benchmark_3d.ugx
GRID_1D=/home/tug41634/two_regions/synapse_simulation_benchmark_1d_soma.ugx
# output
VTK_FOLDER=hybrid_test_par_gmg

# simulation settings and aprameters
SOLVER=GMG # solver id
DT_1D=1e-5 # delta t 1d
DT_3D=1e-4 # delta t 3d
SETTING=all # IP3 and Ryr plus Leak
END_TIME=1.00  # end time
NSYN=100 # 100 synapses

# execution parameters
BINARY=ugshell
NP=56

# execute simulation
mpirun -np "$NP" "$BINARY" \
       -ex "$SCRIPT"  \
       -grid1d "$GRID_1D"  \
       -grid3d "$GRID_3D"  \
       -endTime "$END_TIME" \
       -outName "$VTK_FOLDER" \
       -dt1d "$DT_1D" \
       -dt3d "$DT_3D" \
       -setting "$SETTING" \
       -solver "$SOLVER"  \
       -nSyn "$NSYN" \
       -verbose1d -verbose3d
