#!/bin/bash

set -euo pipefail

# Keep each MPI rank to one CPU thread
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND="cores"         # harmless on a VM
export MPICH_GPU_SUPPORT_ENABLED=0    # harmless if not MPICH
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/diiid-helicon/case_196154/DIII-D_helicon_data

# Show GPU list (optional)
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L || true

# Make only GPU 0 visible (change to 1/2/3 if you want a different GPU)
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=0

mpirun --bind-to core -np 4 "$bin" \
  "${td}/Model5_b-nomesh.osh" \
  "${td}/gitrm_new_4.ptn" \
  "${td}/profilesDIIID.nc" \
  "${td}/surf_fields.txt" \
  "${td}/bfield_diiid.nc" \
  "${td}/lu.nc" \
  "${td}/ADAS_Rates_C.nc" \
  "${td}/ftridynSelf.nc" \
  "${td}/surface_model_GITRm_rustbca_C_W_d.nc" \
  -
