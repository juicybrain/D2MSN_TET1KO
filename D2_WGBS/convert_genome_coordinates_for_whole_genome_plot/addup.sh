#!/usr/bin/bash


add_up(){

   cat $4 |  awk -v OFS="\t" -v chr=${1} -v acL=${3} '{if($1==chr) print $1,$2+ acL,$3+acL, $4}' >>${4}\.new.txt

}



cat ~/ref/genome/male.mm10/male.mm10.genomeSize |grep -v M | awk -v OFS="\t" 'BEGIN{a=0}{a=a+$2}{print $1,$2,a-195471971}' > size_switch.txt


cat /home/yli/ref/genome/chr.addup.txt | while read line
do

    add_up ${line} $1


done
