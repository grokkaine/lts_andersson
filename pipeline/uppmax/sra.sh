#!/bin/bash -l
#SBATCH -A b2014214
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH --mail-user=sergiu.netotea@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -J sample_download
#SBATCH --output=sra%j.out
#SBATCH --error=sra%j.err

module load bioinfo-tools
module load sratools/2.8.0

cd {run_dir}
mkdir qc2
mkdir ./qc2/{sample_name}
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/{sample_name:.3}/{sample_name:.6}/{sample_name}/{sample_name}.sra -O {sample_name}.sra
fastq-dump --split-files ./{sample_name}.sra
/home/sergiu/programs/FastQC/fastqc -o ./qc/{sample_name}/ --nogroup ./{sample_name}_1.fastq ./{sample_name}_2.fastq
cutadapt -m 50 --length-tag "length=" -q 20 -a AGATCGGAAGAG -A AGATCGGAAGAG -o {sample_name}_1.out.fastq -p {sample_name}_2.out.fastq {sample_name}_1.fastq {sample_name}_2.fastq
/home/sergiu/programs/FastQC/fastqc -o ./qc2/{sample_name}/ --nogroup ./{sample_name}_1.out.fastq ./{sample_name}_2.out.fastq
/home/sergiu/programs/fastQValidator/bin/fastQValidator ./{sample_name}_1.out.fastq >> ./qc2/{sample_name}/{sample_name}_1.validator.txt
/home/sergiu/programs/fastQValidator/bin/fastQValidator ./{sample_name}_2.out.fastq >> ./qc2/{sample_name}/{sample_name}_2.validator.txt
