#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/ITER/ITER_antenna_data

mpirun --bind-to core -np 1 $bin \
  ${td}/ITER_geom_antenna.osh \
  ${td}/ITER_antenna_1.ptn \
  ${td}/profiles_iter.nc \
  ${td}/surf_fields.txt \
  ${td}/bfield_iter.nc \
  ${td}/lu.nc \
  ${td}/ADAS_Rates_W.nc \
  ${td}/ftridynSelf.nc \
  ${td}/surface_model_GITRm_rustbca_C_W_d.nc \
  -

