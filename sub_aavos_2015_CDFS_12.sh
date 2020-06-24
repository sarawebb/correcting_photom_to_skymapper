#!/bin/bash 
#SBATCH --job-name=DWF_AAVOS
#SBATCH --output=3hr_AAVOS.out
#SBATCH --error=AAVOS.err
#SBATCH --partition=skylake
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1
#SBATCH --time=128:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --tmp=20GB

python correct_AAVSO_CDFS_2015_12.py
