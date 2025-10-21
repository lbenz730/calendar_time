#!/usr/bin/bash

#SBATCH -c 64 ## number of cores
#SBATCH -t 0-16:00 ## amount of time in D-HH:MM
#SBATCH -p fasse_bigmem ## Partition to submit to
#SBATCH --mem=400000 ## memory pool for all cores
#SBATCH -o logs/EIF_transport/log.stdout_%a ## STDOUT
#SBATCH -e logs/EIF_transport/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-84

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/tv_effects

Rscript scripts/analysis/standardization_EIF_weight.R  $SLURM_ARRAY_TASK_ID
