# -*- coding: utf-8 -*-
"""
Created on Sat May 27 15:12:06 2017

@author: dongp
modified by YLi

"""
import sys
import subprocess
import gc
import time
import os
import gzip

in_name=sys.argv[1]
out_name=sys.argv[2]
suffix=sys.argv[3]
tag=sys.argv[4]
ref='/home/yli/ref/genome/hgT2T/hgT2T'
thread=str(4)
fastq_path='./'
bam_path='./bam/'
tmp_path='./tmp/'
MIN_MAPQ = 11
cst=100
flen=150
step=30

def sleep():
    for _ in range(3):
        time.sleep(0.1)
    gc.collect()
    for _ in range(3):
        time.sleep(0.1) 

def filter_fastq(mapped,fastq0,fastq1):
    with gzip.open(fastq0, 'rt') as infile, gzip.open(fastq1, 'wt') as outfile:
        while True:
            f_id = infile.readline().strip()
            if not f_id:
                break
            
            f_id = f_id[1:]
            f_seq = infile.readline().strip()
            infile.readline()
            f_q = infile.readline().strip()
            read_id = f_id.split()[0]
            
            if read_id not in mapped:
                outfile.write(f"@{f_id}\n{f_seq}\n+\n{f_q}\n")


def iterative_mapping(side):
    fastq0=fastq_path+in_name+side+suffix
    log1=out_name+'_'+side+tag+'.log'
    bam_out=bam_path+out_name+'_'+side+tag
    b1=[]
#    b1.append(subprocess.Popen(['awk',"""{OFS="\\t"; if (1==1) print}"""],
#                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
#                         bufsize=-1))
#    b1.append(subprocess.Popen(['samtools', 'view', '-bS','-'],
#                               stdin=b1[0].stdout,stdout=open(bam_out,'w'),
#                               bufsize=-1,shell=False))
    b1.append(subprocess.Popen(['samtools', 'view', '-bS','-'],stdin=subprocess.PIPE,stdout=open(bam_out,'w'),bufsize=-1,shell=False))

    for i,j in enumerate(range(cst,-1,-step)):
        if i==0:
            mapped=set()
            p1 = subprocess.Popen(['bowtie2','-x',ref, '-U',fastq0,'-3',str(j),
            '-p',thread,'--very-sensitive'],stdout=subprocess.PIPE,stderr=open(log1,'w'))
            for line in p1.stdout:
                line = line.decode('utf-8')
                if line.startswith('@'):
                    b1[0].stdin.write(line.encode('utf-8'))
                    continue
                items=line.split()
                if int(items[4])<MIN_MAPQ:
                    continue
                b1[0].stdin.write(line.encode('utf-8'))
                mapped.add(items[0])
        else:
            sleep()
            fastq1=tmp_path+in_name+str(i)+'_'+str(side)+suffix
            filter_fastq(mapped,fastq0,fastq1)
            mapped=set()
            fastq0=fastq1
            p1 = subprocess.Popen(['bowtie2','-x',ref,'-U',fastq0,'-3',str(j),
                '-p',thread,'--very-sensitive'], stdout=subprocess.PIPE,stderr=open(log1,'a'))
            for line in p1.stdout:
                line = line.decode('utf-8')
                if line.startswith('@'):
                    continue
                items=line.split()
                if int(items[4])<MIN_MAPQ:
                    continue
                b1[0].stdin.write(line.encode('utf-8'))
                mapped.add(items[0])
    b1[0].communicate()
    b1[-1].communicate()
    mapped=set()
    sleep()
#    for i,j in enumerate(range(cst,-1,-step)):
#        if i>0:
#            fastq1=tmp_path+in_name+str(i)+'_'+str(side)+suffix
#            os.remove(fastq1)
            
for n in range(2):
    iterative_mapping(str(n+1))    

