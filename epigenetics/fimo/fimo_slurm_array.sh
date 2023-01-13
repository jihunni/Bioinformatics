#!/bin/bash
#
#SBATCH --job-name=fimo
#
#SBATCH --account=jihun
#SBATCH --mail-user=jihun@gm.gist.ac.kr
#SBATCH --mail-type=BEGIN,FAIL,END
#
#
#SBATCH --ntasks=1
#SBATCH --array=0-60
#SBATCH --cpus-per-task=1
#SBATCH --mem 4G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_result_%j.txt
#SBATCH -o %x_%A.o%j
#SBATCH -e %x_%A.e%j

# stdout : basic info
echo "SLURM running info"
echo "Start Time        = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

module purge
module load MEME_Suite/5.4.1

~/data/motif/slurm_array/fimo.sh
sleep 1 # pause to be kind to the scheduler

echo "Finish Time       = $(date)"
