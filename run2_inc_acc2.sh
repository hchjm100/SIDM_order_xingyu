#!/bin/bash -l
#
#SBATCH --output=Gadget2SIDM_run2_inc_acc2.out     # Standard output and error log
#
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=128              # Total MPI tasks (matching mpirun -np #)
#SBATCH --cpus-per-task=1       # Number of CPU cores per MPI task
#SBATCH --mem=200G                # Total memory per node (adjust as needed)
#SBATCH --time=30-00:00:00       # Time limit (D-HH:MM:SS) this case: 10 days
#SBATCH --mail-user=hsang012@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="acc_inc2"
#SBATCH -p epyc                 # Partition (change to epyc, intel, batch, etc.)

module load fftw/2.1.5

# 2) Print some diagnostic info
echo "Starting Gadget2SIDM job on ${HOSTNAME}"
echo "Date: $(date)"
echo "Directory: $(pwd)"

# 3) Run your program
# If the executable and .param file are in the current directory:
mpirun -np 128 ./GadCrater3 run2_inc_acc2.param

# 4) Print completion time
echo "Finished at: $(date)"

