#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/WEST/Guillaume_runs/1p1e17/WEST_data

mpirun --bind-to core -np 1 $bin \
  ${td}/Model2-case2.osh \
  ${td}/WEST_1.ptn \
  ${td}/profilesWEST_ne_lim_7p3e18.nc \
  ${td}/surf_fields_ne_lim_7p3e18.txt \
  ${td}/BfieldWEST.nc \
  ${td}/lu.nc \
  ${td}/ADAS_Rates_W.nc \
  ${td}/ftridynSelf.nc \
  ${td}/surface_model_GITRm_rustbca_C_W_d.nc \
  -

