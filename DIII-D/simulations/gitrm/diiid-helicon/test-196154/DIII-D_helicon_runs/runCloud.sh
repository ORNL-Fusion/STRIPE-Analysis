#!/bin/bash

export SLURM_CPU_BIND="cores"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/diiid-helicon/case_196154/DIII-D_helicon_data

mpirun --bind-to core -np 1 $bin \
  ${td}/Model5_b-nomesh.osh\
  ${td}/DIII-D_helicon_1.ptn \
  ${td}/profilesDIIID.nc \
  ${td}/surf_fields.txt \
  ${td}/bfield_diiid.nc \
  ${td}/lu.nc \
  ${td}/ADAS_Rates_C.nc \
  ${td}/ftridynSelf.nc \
  ${td}/surface_model_GITRm_rustbca_C_W_d.nc \
  -

