#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/WEST/WEST_data

mpirun --bind-to core -np 1 $bin \
  ${td}/Model2-case2.osh \
  ${td}/WEST_1.ptn \
  ${td}/profilesWEST.nc \
  ${td}/surf_points_rot.txt \
  ${td}/BfieldWEST.nc \
  ${td}/lu.nc \
  ${td}/ADAS_Rates_W.nc \
  ${td}/ftridynSelf.nc \
  ${td}/surface_model_GITRm_rustbca_C_W_d.nc \
  -

