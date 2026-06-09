#!/bin/bash
#SBATCH --job-name=ieads_west_w/PetraM
#SBATCH --output=ieads.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=168:05:00
#SBATCH -p regular
#SBATCH --mail-user=kumara@ornl.gov
#SBATCH --mail-type=ALL
# ----------------------------
#srun python runGITR.py
srun  python3 runGITR.py
