#!/bin/bash
#needs 1 arguments
#$1 path to demultiplexed files
#$2 Ext of the raw files fastq or fastq.gz
#$3 FW primer (for read $sample_1.fastq.gz)
#$4 Require MINLENGTH overlap between read and adapter for an adapter to be found
#$5 extra params for cutadapt
mkdir -p  $1/primer_removed 
for fw in $1/*.$2 ; do
   echo $sample
   ext=$2
   sample=$(echo $fw | awk -F"/" '{print $NF}' | sed "s/.$ext//");
   cutadapt -g "$3"  -O $4 $5 -o $1/primer_removed/${sample}_1..fastq.gz $fw  > $1/$sample.cutadapt.log
   res=$(grep "(passing filters)" $1/primer_removed/$sample.cutadapt.log | awk '{print $5"\t"$6}')
   echo -e $sample"\t"$res >> $1/summary.txt
done

