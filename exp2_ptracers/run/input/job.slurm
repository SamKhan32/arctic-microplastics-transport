#!/bin/bash
#SBATCH --job-name=arctic_circ
#SBATCH --time=048:00:00
#SBATCH --ntasks=90
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sphill86@charlotte.edu
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

echo "======================================================"
echo "Start Time : $(date)"
echo "Submit Dir : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Partition : $SLURM_JOB_PARTITION"
echo "Num Tasks : $SLURM_NTASKS total"
echo "======================================================"

module purge
module load gcc/13.2.0
module load openmpi/4.1.6

set -euo pipefail
cd /users/sphill86/baseline_exp_llc270/run/execution
export OMP_NUM_THREADS=1

echo "Launching MITgcm with $SLURM_NTASKS MPI ranks..."
srun -n 90 ./mitgcmuv
cd /users/sphill86/baseline_exp_llc270/run/execution
tar -cf ~/arctic_output_2years.tar *.data *.meta

echo "======================================================"
echo "End Time : $(date)"
echo "======================================================"
