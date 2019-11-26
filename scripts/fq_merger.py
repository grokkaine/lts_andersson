"""
cd /proj/b2010008/nobackup/projects/crispr/sergiu/src/ && \
rm fq_merger.py && nano fq_merger.py

module use /proj/b2013006/sw/modules && \
module load miniconda3 && \
source activate andersson


"""

import sys

def merge_chunks(tmerge, nmerge, merge10, samples):
    rem_samples = set(["P1994_105_R1","P1994_105_R2"])
    from Bio import SeqIO
    import gzip
    handle_merge10 = gzip.open(merge10, "wt")
    # First add the reads from the new merged file to the surface sample merge
    handle_nmerge = gzip.open(nmerge, "rt")
    fq = SeqIO.parse(handle_nmerge, "fastq")
    for read in fq:
        sample = read.id.split(":")[0]
        if sample in samples:
            handle_merge10.write(read.format("fastq"))
    handle_nmerge.close()
    # Now add the remaining reads from the old merged file and add the old to the new
    handle_tmerge = gzip.open(tmerge, "rt")
    handle_nmerge = gzip.open(nmerge, "at")
    fq = SeqIO.parse(handle_tmerge, "fastq")
    for read in fq:
        sample = read.id.split(":")[0]
        if sample in samples:
            handle_merge10.write(read.format("fastq"))
        if not sample in rem_samples:
            handle_nmerge.write(read.format("fastq"))
    handle_tmerge.close()
    handle_nmerge.close()
    handle_merge10.close()
    return


def run():
    samples10 = ["P1994_101", "P1994_104", "P1994_107", "P1994_110", "P1994_113", "P1994_116", "P1994_119", "P1994_122", "P1994_125", "P1994_128"]
    samples = set([s+"_R1" for s in samples10] + [s+"_R2" for s in samples10])
    data_dir = "/proj/b2010008/nobackup/projects/crispr/data/transect/"
    tmerge = data_dir + "temp_transect_merged.fastq.gz"
    nmerge = data_dir + "transect_merged.fastq.gz"
    merge10 = data_dir + "transect10_merged.fastq.gz"
    merge_chunks(tmerge, nmerge, merge10, samples)
    print("Done!")
    return

run()
