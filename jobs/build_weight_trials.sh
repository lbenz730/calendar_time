#!/usr/bin/bash

#SBATCH -c 4 ## number of cores
#SBATCH -t 0-00:30 ## amount of time in D-HH:MM
#SBATCH -p serial_requeue ## Partition to submit to
#SBATCH --mem=170000 ## memory pool for all cores
#SBATCH -o logs/build_weight/log.stdout_%a ## STDOUT
#SBATCH -e logs/build_weight/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-84

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/tv_effects

Rscript scripts/data/build_weight_trials.R  $SLURM_ARRAY_TASK_ID
