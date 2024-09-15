### Processing and alignment of Hi-C sequencing reads
#
#  Author Yuxiang Li
#

# Trim adaptor
trim_galore --paired   ${spl}\_R1.fastq.gz ${spl}\_R2.fastq.gz

#  Align reads with Hi-C Pro with the following parameters in config-hicpro.txt

	nohup HiC-Pro -i  /home/yli/workspace/YuxiangTempData/HiC/D2_HiC/rawdata -o ./D2_HiC -c config-hicpro.arima.txt &

# Create ".hic" files for JuicerTools HiC-Pro valid pairs 
	   for VALIDPAIRS in `ls *Pairs`; do  awk '{$4=$4!="+"; $7=$7!="+"; n1=split($9, frag1, "_"); n2=split($10, frag2, "_"); } $2<=$5{print $1, $4, $2, $3, frag1[n1], $7, $5, $6, frag2[n2], $11, $12 }$5<$2{ print $1, $7, $5, $6, frag2[n2], $4, $2, $3, frag1[n1], $12, $11}' $VALIDPAIRS | LANG=C sort -T tmp -k3,3d -k7,7d -S 50% --parallel=4 > tmp/${VALIDPAIRS}\_allValidPairs.pre_juicebox_sorted & done
     	nohup java -Xmx128g -jar ~/src/juicertools/juicer_tools_1.22.01.jar pre  -j 32  tmp/JMC.allValidPairs_allValidPairs.pre_juicebox_sorted  ./test_2/JMC.D2.noK.hic   ~/ref/genome/male.mm10/male.mm10.natural.genomeSize  &
     	nohup java -Xmx128g -jar ~/src/juicertools/juicer_tools_1.22.01.jar pre  -j 32  tmp/YL.allValidPairs_allValidPairs.pre_juicebox_sorted  ./test_2/YL.D2.noK.hic   ~/ref/genome/male.mm10/male.mm10.natural.genomeSize  &

	

# identify TAD at 25k resolution for merged four groups Hi-C data: male Con, male KO, female Con, female KO
java -jar  juicer_tools_1.22.01.jar arrowhead -r 25000 -k KR --ignore-sparsity ${spl}\.hic ${spl}\.25k.KR

# Identify loop at 10000, 25000, 50000 resolution
java -jar juicer_tools_1.22.01.jar hiccups  -k KR -r 5000,10000,25000 -f 0.1 --ignore-sparsity ${spl}  ${spl}\_10.25.50k --cpu
java -jar juicer_tools_1.22.01.jar hiccups  -k KR -r 50000 -f 0.1 --ignore-sparsity ${spl}  ${spl}\_10.25.50k --cpu

# create Homer tag depository
  cut -f 1-7 validPairs > validPairs.homer
  makeTagDirectory ${name} -format HiCsummary ${spl}
# calculate the loop score for each Hi-C replicate
findTADsAndLoops.pl score  -loop ../D1_loop.bedpe  -o male.loop   -d <Hi-C replicates directory> -cpu 4 -res 25000  -window 50000 -raw -normTotal 1e8 



# convert ".hic" file to ".cool" file using Hicexplorer
hicConvertFormat   --matrices ${spl}  --outFileName ${spl}\.mcool --inputFormat hic --outputFormat cool

# get matrix of certain resolution
hicConvertFormat --matrices ${spl}\.mcool::/resolutions/25000 --outFileName ${spl}\_250k_${i}\.cool --inputFormat cool --outputFormat cool 

# mask bad regions (downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz and http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/genomicSuperDups.txt.gz )
hicAdjustMatrix -m ${spl} --action mask --regions   male.mm10.badrigion.GT25k.bed  -o mask/${spl}\.mask.cool

# normalize the reads of Hi-C replicates
current_time=$(date +%Y.%m.%d-%H.%M.%S)
hicNormalize -m  <all replicates>  --normalize smallest  -o <normalized files>  2>${current_time}\.normalTosmallest.log  

# Compare Tet1 KO to Tet1 Con Hi-C matrices
hicCompareMatrices -m ${spl}\_250k.cool.mask.nrl.cool F_N_250k.cool.mask.nrl.cool --outFileName ${spl}\_nrl_PosMinusNeg.log2ratio.cool --operation log2ratio

# KR correction of Hi-C matrices
hicCorrectMatrix correct --matrix  ${spl}\_25k.mask.cool   --correctionMethod  KR -o ${spl}\_25k.mask.KR.cool

# Hi-C A/B compartment analysis PC1 analysis with hicPCA, K27ac ChIPseq data as marker for active chromatin

hicPCA --matrix ../${spl}\.2Rep.allValidPairs.hic.mcool::/resolutions/1000000    --outputFileName  ${spl}\.distalN.1M.pca.txt  --whichEigenvectors 1  --format bedgraph --method dist_norm  --extraTrack  GSE197668_MSD1K27A.bw  --histonMarkType active
hicPCA --matrix ../${spl}\.2Rep.allValidPairs.hic.mcool::/resolutions/1000000    --outputFileName  ${spl}\.distalN.1M.pca.txt  --whichEigenvectors 1  --format bedgraph --method dist_norm  --extraTrack  GSE197668_FSD1K27A.bw  --histonMarkType active







