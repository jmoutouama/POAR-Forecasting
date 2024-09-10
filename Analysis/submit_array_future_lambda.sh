#!/bin/bash
 
#SBATCH --job-name=test_job
#SBATCH --chdir=/work/compagna
#SBATCH --output=/work/%u/%x-%A-%a.out
#SBATCH --time=0-00:60:00
#SBATCH --mem-per-cpu=8G

module load GCC/12.2.0 
module load OpenMPI/4.1.4
module load R/4.2.2
module spider renv

CLIMATE="$2"

app \
  --input=/data/project/input/$SLURM_ARRAY_TASK_ID \
  --output=/work/$USER/$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID
  
Rscript --vanilla /home/$USER/poar/future_clust_lambda.R \
	"$SLURM_ARRAY_TASK_ID" \
	"$CLIMATE"