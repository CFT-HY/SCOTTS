#!/bin/bash -l
#SBATCH -J gc-test
#SBATCH -e gc-test_err_%j
#SBATCH -o gc-test_out_%j
#SBATCH --mem-per-cpu=2000
#SBATCH -t 36:00:00
#SBATCH -n 256
#SBATCH -p parallel

module load intel

mkdir -p /wrk/weir/gc-hydro/gc-test/

srun ~/appl_taito/projects/gc-hydro/run/hydro input-intermediate-0.4-both
