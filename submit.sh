#!/bin/bash -l

#SBATCH -J SIDM50-iso-m3e8-c28.4-low
#SBATCH -o nohup.out
#SBATCH -p epyc
#SBATCH --ntasks=128
#SBATCH --nodes=1
#SBATCH --mem=200gb
#SBATCH --time=168:00:00

# Load needed modules
# You could also load frequently used modules from within your ~/.bashrc

# Run job utilizing all requested processors
# Please visit the namd site for usage details: http://www.ks.uiuc.edu/Research/namd/

module load fftw/2.1.5

mpirun -np 128 ./GadCrater3 StelPert.param 


