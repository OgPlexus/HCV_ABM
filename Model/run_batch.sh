#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH -t 4:00:00
#SBATCH --array=1-30
#SBATCH --constraint=skylake
#SBATCH --account=ccmb-condo

module load python/3.5.2

python3 HCV_Sim_v2_0.py $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID

sbatch --output=/dev/null --error=/dev/null
