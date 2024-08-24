# -*- coding: utf-8 -*-
"""
2024.08.24

Yuxiang Li

"""

import subprocess
import sys
import os
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def split_id(seq_id):
    return seq_id.split()[0]

# Get the input file handle from command line arguments
fh = sys.argv[1]

# Open output files for writing compressed data using pigz
ofh1 = open(fh + '.clean.R1.fq.gz', 'wb')  # Use binary mode 'wb' for writing gzip files
ofh2 = open(fh + '.clean.R2.fq.gz', 'wb')  # Use binary mode 'wb' for writing gzip files
out1 = subprocess.Popen(["pigz", "-c"], stdin=subprocess.PIPE, stdout=ofh1, bufsize=-1)
out2 = subprocess.Popen(["pigz", "-c"], stdin=subprocess.PIPE, stdout=ofh2, bufsize=-1)

# Open the gzip FASTQ files for reading using gzip.open
f_fastq = FastqGeneralIterator(gzip.open(fh + '.R1.fq.gz', 'rt'))
r_fastq = FastqGeneralIterator(gzip.open(fh + '.R2.fq.gz', 'rt'))

n = 0

for (f_id, f_seq, f_q) in f_fastq:
    match_id_f = split_id(f_id)

    # Reset the r_fastq iterator to start from the beginning for each query
    r_fastq = FastqGeneralIterator(gzip.open(fh + '.R2.fq.gz', 'rt'))  # Re-open r_fastq for each f_fastq read

    for (r_id, r_seq, r_q) in r_fastq:
        match_id_r = split_id(r_id)

        # Check if both reads are paired by matching their IDs
        if match_id_f == match_id_r:
            out1.stdin.write(f"@{f_id}\n{f_seq}\n+\n{f_q}\n".encode('utf-8'))
            out2.stdin.write(f"@{r_id}\n{r_seq}\n+\n{r_q}\n".encode('utf-8'))
            n += 1
            break  # Stop searching in r_fastq once a match is found

print(fh + '\t' + str(n) + '\n')

out1.communicate()
out2.communicate()
