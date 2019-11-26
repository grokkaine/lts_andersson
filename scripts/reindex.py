from Bio import SeqIO

import sys

def test(stdin, stdout, add="_addthis"):
    fo = open(stdout, 'wt')
    with open(stdin, 'rt') as f:
        for l in f:
            print("reading line:")
            fo.write(l+add)
    return


def reindex_fq(input_stream, output_stream, sid):
    infile = open(input_stream, "rU")
    outfile = open(output_stream, "w")
    for read in SeqIO.parse(input_stream, "fastq"):
        read.id = sid + ":" + read.id
        outfile.write(read.format("fastq"))
    infile.close()
    outfile.close()
    return


reindex_fq(sys.argv[1],sys.argv[2],sys.argv[3])
