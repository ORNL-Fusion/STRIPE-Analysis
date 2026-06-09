#!/bin/bash
#SBATCH -A m77_g
#SBATCH -C gpu
#SBATCH -G 1
##SBATCH -J ieads
#SBATCH -q regular
#SBATCH --mail-user=kumara@ornl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 128
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=none

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load cmake
module load cray-hdf5
module load cray-netcdf
module load python


#srun python runGITR.py
srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=single:1 python runGITR.py
