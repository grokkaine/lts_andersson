#!/bin/bash -l
#SBATCH -A b2014214
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --mail-user=sergiu.netotea@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -J sample_merger
#SBATCH --output=sample_merger-%j.out
#SBATCH --error=sample_merger-%j.err

# use - p devcore, -t 00:05:00 when testing
# use -p core, -t 3-00:00:00 when not testing

# use the crispr virtual environment
source /home/sergiun/crispr/bin/activate
python ../pipeline.py /proj/b2010008/nobackup/projects/crispr/sergiu/src/configure.yaml
# get out of the virtual environment
deactivate
