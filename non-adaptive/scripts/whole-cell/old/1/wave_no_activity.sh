#!/bin/sh
#PBS -l walltime=06:00:00
#PBS -V
#PBS -m bea
#PBS -N wave3d_revised_multi_measure_zones_NO_ACTIVITY
#PBS -q normal
#PBS -l nodes=1:ppn=1
cd "$PBS_O_WORKDIR"
SCRIPT=wave3d_revised_multi_measure_zones_NO_ACTIVITY.lua
GRID=/home/tug41634/single_branch2_w_branch_meas_new.ugx
SOLVER=GMG
mpirun -np 1  ugshell -ex "$SCRIPT" -grid "$GRID" -dt 1e-6 -endTime 100.00 -solver "$SOLVER" -tol 0.02 -ryrDensity 10 -setting ryr -vtk -pstep 0.0 -outName /home/tug41634 -vtkFolder vtk8 -suffix vtk8/meas/ -measSuffix meas -vtk
