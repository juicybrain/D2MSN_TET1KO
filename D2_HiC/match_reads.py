"""
0824.2024
Yuxiang Li
"""
import gzip
import subprocess
import sys,os
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def split_id(seq_id):
    return seq_id.split()[0]


fh=sys.argv[1]

ofh1=open(fh+'.clean.R1.fq.gz','wb')
ofh2=open(fh+'.clean.R2.fq.gz','wb')
out1=subprocess.Popen(["pigz", "-c"], stdin = subprocess.PIPE, \
stdout = ofh1, bufsize = -1)
out2=subprocess.Popen(["pigz", "-c"], stdin = subprocess.PIPE, \
stdout = ofh2, bufsize = -1)
#f_fastq = FastqGeneralIterator(open(fh+'.R1.fq.gz'))
#r_fastq = FastqGeneralIterator(open(fh+'.R2.fq.gz'))

f_fastq = FastqGeneralIterator(open(fh + '.R1.fq', 'rt'))  # Use 'rt' mode for reading text from gzip
r_fastq = FastqGeneralIterator(open(fh + '.R2.fq', 'rt'))  # Use 'rt' mode for reading text from gzip


r_index = SeqIO.index_db('tmp.idx', fh+'.R2.fq', "fastq")
n=0
for (f_id, f_seq, f_q) in f_fastq:
    match_id = split_id(f_id)
    if match_id in r_index:
        r_read = r_index.get_raw(match_id)
        out1.stdin.write(f"@{f_id}\n{f_seq}\n+\n{f_q}\n".encode('utf-8'))   
        out2.stdin.write(r_read)
        n+=1
print(fh+'\t'+str(n)+'\n')
out1.communicate()
out2.communicate()
os.remove('tmp.idx')
os.remove(fh+'clean.R1.fastq')
os.remove(fh+'clean.R2.fastq')
