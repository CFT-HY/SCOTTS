#!/bin/bash -l
#SBATCH -J rerun-0.4
#SBATCH -e rerun_err_%j
#SBATCH -o rerun_out_%j
#SBATCH -t 5:00:00
#SBATCH -N 30
#SBATCH -p parallel
#SBATCH --no-requeue

module load PrgEnv-intel

(( ncores = $SLURM_NNODES * 24))

mkdir -p /wrk/weir/rerun-0.4/

aprun -n $ncores -N 24 ../hydro input-intermediate-0.4-both
