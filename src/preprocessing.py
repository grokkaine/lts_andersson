"""
This file contains different file processing commands.
- merge all FASTQ files in a given list for use by Crass
- download the A SRA archives
"""


def run_batch(location, template_fpath, run_fpath, params):
    """
    The parameters are given in a name: value dictionary and are inserted in
    the template script which is then submitted to UPPMAX.
    """
    import os
    f = open(template_fpath, 'r')
    script = "".join(f.readlines())
    f.close()
    with open(run_fpath, 'w') as f:
        f.write(script.format(**params))
    if location == "UPPMAX":
        os.system("sbatch " + run_fpath)
    elif location == "workpc":
        # import os
        # import stat
        # st = os.stat('somefile')
        # os.chmod('somefile', st.st_mode | stat.S_IEXEC)
        import subprocess
        subprocess.call(["chmod", "a+x", run_fpath])
        subprocess.call([run_fpath])
    return


def merge_fq(filename, file_list, sample_list):
    import gzip
    import datetime
    fo = gzip.open(filename, 'wt')
    with open("merge_log.txt", 'a') as f:
        f.write('Timestamp: {:%Y-%m-%d %H:%M:%S}'
                .format(datetime.datetime.now()) + "\n")
        f.write("Output to" + filename + "\n")
    for samplefilen, sid in zip(file_list, sample_list):
        with open("merge_log.txt", 'a') as f:
            f.write('Timestamp: {:%Y-%m-%d %H:%M:%S}'
                    .format(datetime.datetime.now()) + "\n")
            f.write("Adding file:" + samplefilen + "\n")
        with gzip.open(samplefilen, 'rt') as fin:
            for line in fin:
                if line[0] == '@':
                    line = '@' + sid + ":" + line[1:]
                fo.write(line)
    fo.close()
    return


def reindex_fq(input_stream, output_stream, sid):
    from Bio import SeqIO
    fq = SeqIO.parse(input_stream, "fastq")
    for read in fq:
        read.id = sid + ":" + read.id
        output_stream.write(read.format("fastq"))
    return


def test(filen):
    with open(filen, 'w') as f:
        f.write("test")
    print("Done test!")
    return


def SRA_download(template_fpath, run_dir, download_loc, sample_names,
                 location):
    """
    The ftp download uses the documentation available here:
    https://www.ncbi.nlm.nih.gov/books/NBK158899/
    I download the samples and unpack them straight on UPPMAX.
    """
    for sample_name in sample_names:
        run_batch(location, template_fpath, run_dir+'/'+sample_name+'.sh',
                  {'sample_name': sample_name, 'run_dir': run_dir})
    print("Submitted the SRA download batch!")
    return
