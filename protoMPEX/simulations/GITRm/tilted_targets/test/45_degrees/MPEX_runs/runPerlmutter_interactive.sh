#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/global/homes/a/atul19/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/pscratch/sd/a/atul19/gitrRuns/MPEX/tilted_targets/45_degrees/MPEX_data

srun $bin \
  ${td}/Model_45.osh \
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
