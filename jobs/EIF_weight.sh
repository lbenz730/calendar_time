#!/usr/bin/bash

#SBATCH -c 64 ## number of cores
#SBATCH -t 0-24:00 ## amount of time in D-HH:MM
#SBATCH -p fasse_bigmem ## Partition to submit to
#SBATCH --mem=190000 ## memory pool for all cores
#SBATCH -o logs/EIF_weight/log.stdout ## STDOUT
#SBATCH -e logs/EIF_weight/log.stderr ## STDERR
#SBATCH --account=haneuse_lab

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/tv_effects

Rscript scripts/analysis/EIF_weight_analysis.R
Rscript scripts/analysis/EIF_weight_projections.R
