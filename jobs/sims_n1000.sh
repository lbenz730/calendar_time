#!/usr/bin/bash

#SBATCH -c 16 ## number of cores
#SBATCH -t 0-01:00 ## amount of time in D-HH:MM
#SBATCH -p serial_requeue ## Partition to submit to
#SBATCH --mem=30000 ## memory pool for all cores
#SBATCH -o logs/sims_n1000/log.stdout_%a ## STDOUT
#SBATCH -e logs/sims_n1000/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-1000

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/tv_effects

Rscript scripts/simulations/run_simulation_pipeline.R  $SLURM_ARRAY_TASK_ID $* 1000
