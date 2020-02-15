#!/bin/bash
#SBATCH --account=pi-larogers
#SBATCH --job-name=MalskyTest
#SBATCH --time=36:00:00
##SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=isaacmalsky@gmail.com

module load Anaconda3/5.1.0

export OMP_NUM_THREADS=16

./clean
./clear_inlists
./mk
python run_orion.py

