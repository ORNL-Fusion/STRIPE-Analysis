#!/bin/bash

set -euo pipefail

# Keep each MPI rank to one CPU thread
export OMP_NUM_THREADS=1
export SLURM_CPU_BIND="cores"         # harmless on a VM
export MPICH_GPU_SUPPORT_ENABLED=0    # harmless if not MPICH
export MPICH_ABORT_ON_ERROR=1
ulimit -c unlimited

bin=/home/cloud/myRepos/gitrm_root/build-gitrm-perlmutter-cuda/GITRm
td=/home/cloud/data/GITRm/diiid-helicon/case_200882/DIII-D_helicon_data

# Show GPU list (optional)
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L || true

# Make per-rank GPU selection (rank->GPU: 0,1,2,3); keep device order stable
export CUDA_DEVICE_ORDER=PCI_BUS_ID
# DO NOT set CUDA_VISIBLE_DEVICES globally anymore
# export CUDA_VISIBLE_DEVICES=0   # <-- removed

# minimal wrapper to map local rank -> GPU ID
wrapper="$(mktemp)"
cat > "$wrapper" <<'EOF'
#!/bin/bash
lr=${MPI_LOCALRANKID:-${OMPI_COMM_WORLD_LOCAL_RANK:-0}}
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=${lr}
exec "$@"
EOF
chmod +x "$wrapper"

mpirun --bind-to core -np 3 "$wrapper" "$bin" \
  "${td}/Model6_b_nomesh.osh" \
  "${td}/gitrm_new_3.ptn" \
  "${td}/profilesDIIID_200882.nc" \
  "${td}/surf_fields.txt" \
  "${td}/bfield_diiid.nc" \
  "${td}/lu.nc" \
  "${td}/ADAS_Rates_C.nc" \
  "${td}/ftridynSelf.nc" \
  "${td}/surface_model_GITRm_rustbca_C_W_d.nc" \
  -

rm -f "$wrapper"
