#!/bin/bash
#SBATCH -J gitrm_1gpu
#SBATCH -A m3236
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus=1
#SBATCH --gpu-bind=none
#SBATCH --mail-user=kumara@ornl.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o gitrm_%j.out
#SBATCH -e gitrm_%j.err

set -euo pipefail

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/global/homes/a/atul19/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/pscratch/sd/a/atul19/gitrRuns/tilted_targets/45_degrees/MPEX_data

# ---- minimal automation: robust CPUs-per-task (SLURM var isn't always set) ----
CPT=${SLURM_CPUS_PER_TASK:-16}

# Run on 1 GPU, 1 task
srun --ntasks=1 --cpus-per-task=$CPT --gpus=1 \
  "$bin" \
  "${td}/Model_45_split.osh" \
  "${td}/MPEX_split_1.ptn" \
  "${td}/profilesProtoMPEX.nc" \
  "${td}/surf_fields_New.txt" \
  "${td}/BfieldProtoMPEX.nc" \
  "${td}/lu.nc" \
  "${td}/ADAS_Rates_W.nc" \
  "${td}/ADAS_Rates_W.nc" \
  "${td}/ftridynSelf.nc" \
  "${td}/TaW_surface_response_v2.nc" \
  -
