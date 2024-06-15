#!/bin/bash
#SBATCH --mem=96G       # Maximum memory required per CPU (in megabytes)
#SBATCH --ntasks-per-node=24
#SBATCH --error=/work/teckyongtan/tecktan/Bundling/Data/Temp/Error_Output/err_%A_%a
#SBATCH --output=/work/teckyongtan/tecktan/Bundling/Data/Temp/Error_Output/out_%A_%a

module load anaconda
conda activate R-CUSTOM

cd /work/teckyongtan/tecktan/Bundling/familyenrollment

Rscript estimation/estimation.R $SLURM_ARRAY_TASK_ID $SLURM_NTASKS_PER_NODE 