#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root_linear_test/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/MPEX/tilted_targets/test/0_degrees/MPEX_data

mpirun --bind-to core -np 1 $bin \
  ${td}/Model_00_hang_surf.osh \
  ${td}/MPEX_split_1.ptn \
  ${td}/profilesProtoMPEX.nc \
  ${td}/surf_fields_New.txt \
  ${td}/BfieldProtoMPEX.nc \
  ${td}/lu.nc \
  ${td}/ADAS_Rates_W.nc \
  ${td}/ADAS_Rates_W.nc \
  ${td}/ftridynSelf.nc \
  ${td}/TaW_surface_response_v2.nc \
  -

