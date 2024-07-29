#!/usr/bin/bash

add_up(){
       # Read parameters into named variables for clarity
       chr=$1
              acL=$3
                 inputFile=$4

                    cat "$inputFile" |  awk -v OFS="\t" -v chr="${chr}" -v acL="${acL}" '{if($1==chr) print $1,$2+acL,$3+acL, $4}' >> ${inputFile}.new.txt
}

# Pre-processing step: Generate size_switch.txt based on genome size information, genomesize file have chrs order by 1..19,X, and Y
cat ~/ref/genome/male.mm10/male.mm10.natural.genomeSize | grep -v M | awk -v OFS="\t" 'BEGIN{a=0}{a=a+$2}{print $1,$2,a-$2}' > chr_switch.txt

# Main processing loop: For each .bed file in the current directory
for spl in *.bed
    do
    # Read each line from chr.addup.txt and call the add_up function with the current line and .bed file as arguments
    cat chr_switch.txt | while read line
    do
        add_up ${line} $spl
done
done
