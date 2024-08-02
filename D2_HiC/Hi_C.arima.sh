###
python digest_genome.py -r ^GATC G^AATC G^ATTC G^ACTC G^AGTC -o female.mm10.arima_GATC_GANTC.bed female.mm10.fasta

###
nohup HiC-Pro -i  /home/yli/workspace/YuxiangTempData/HiC/D2_HiC/rawdata -o ./D2_HiC -c config-hicpro.arima.txt &

