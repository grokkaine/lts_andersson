#!/bin/bash
#SBATCH -A b2010008
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J run_crass
#SBATCH -o run_crass.out
#SBATCH -e run_crass.err

module use /proj/b2013006/sw/modules
module load miniconda3
source activate andersson

cd /proj/b2010008/nobackup/projects/crispr/sergiu/bin/bin && \
./crass -o /proj/b2010008/nobackup/projects/crispr/data/asko/crass \
-l 3 /proj/b2010008/nobackup/projects/crispr/data/asko/merged.fastq.gz >run.log 2>&1

source deactivate andersson

