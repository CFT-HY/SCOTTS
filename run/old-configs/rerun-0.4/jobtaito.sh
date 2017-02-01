#!/bin/bash -l
#SBATCH -J rerun-0.4
#SBATCH -e rerun-0.4_err_%j
#SBATCH -o rerun-0.4_out_%j
#SBATCH --mem-per-cpu=2000
#SBATCH -t 0:30:00
#SBATCH -n 64
#SBATCH -p parallel

module load intel

mkdir -p /wrk/weir/gc-hydro/rerun-0.4/

srun ~/appl_taito/projects/gc-hydro/run/hydro input-intermediate-0.4-both
