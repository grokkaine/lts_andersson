#! /bin/bash
#SBATCH -A b2010008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --mail-user=sergiu.netotea@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -J askomain
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# not needed, I installed it in the project conda
#module load bioinfo-tools
#module load FastQC/0.11.5
#module load sratools/2.8.0
#module load cutadapt/1.12
#module load snakemake/3.8.0

module use /proj/b2013006/sw/modules
module load miniconda3

#conda create --name andersson python=3.5 pyyaml
#source activate andersson
#conda install -c bioconda snakemake
#conda install -c bioconda cutadapt
#conda install -c bioconda fastqc=0.11.5
#conda install -c biobuilds sra-tools=2.8.0
#source deactivate andersson

#http://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html#cluster-configuration
#https://wabi-wiki.scilifelab.se/display/KB/Running+snakemakelib+@uppmax
# TODO:Make a test run

source activate andersson
snakemake -s Snakefile -j 99 --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} -t {cluster.time} -p {cluster.partition} -n {cluster.n}"
source deactivate andersson
