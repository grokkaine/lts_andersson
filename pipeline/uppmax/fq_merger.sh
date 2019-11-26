#!/bin/bash -l
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 9-00:00:00
#SBATCH --mail-user=sergiu.netotea@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -J fq_merger
#SBATCH --output=fq_merger-%j.out
#SBATCH --error=fq_merger-%j.err

cd /proj/b2010008/nobackup/projects/crispr/data/transect/ && \
module use /proj/b2013006/sw/modules && \
module load miniconda3 && \
source activate andersson

python3.5 /proj/b2010008/nobackup/projects/crispr/sergiu/src/fq_merger.py

source deactivate andersson
