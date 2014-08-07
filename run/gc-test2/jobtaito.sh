#!/bin/bash -l
#SBATCH -J gc-test2
#SBATCH -e gc-test2_err_%j
#SBATCH -o gc-test2_out_%j
#SBATCH --mem-per-cpu=3800
#SBATCH -t 3-00:00:00
#SBATCH -n 256
#SBATCH -p parallel

module load intel

mkdir -p /wrk/weir/gc-hydro/gc-test2/

srun ~/appl_taito/projects/gc-hydro/run/hydro input-intermediate-0.4-both
