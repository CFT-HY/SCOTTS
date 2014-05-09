#!/bin/bash -l
#SBATCH -J rerun-0.4
#SBATCH -e rerun_err_%j
#SBATCH -o rerun_out_%j
#SBATCH -t 0:30:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -p test

module load intel

mkdir -p /wrk/weir/rerun-0.4/

aprun -n 16 ~/appl_sisu/projects/newhydro/run/hydro input-intermediate-0.4-both
