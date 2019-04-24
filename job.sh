#!/bin/sh
#SBATCH -p devel
#SBATCH -N 2
#SBATCH --tasks-per-node 8
#SBATCH -t 00:10:00
srun ./src/snap_mkl_vml --fi qasnap/sample/inp --fo snap_test.out
