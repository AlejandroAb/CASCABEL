#!/bin/bash
#needs 1 arguments
#$1 path to demultiplexed files
#$2 extension fastq.gz or only fastq
#$3 FW primer (for read $sample_1.fastq.gz)
#$4 RV primer (for read $sample_2.fastq.gz)
#$5 Require MINLENGTH overlap between read and adapter for an adapter to be found
#$6 extra params for cutadapt 
for fw in $1/*_1.$2 ; do
   echo $sample
   ext=$2
   sample=$(echo $fw | awk -F"/" '{print $NF}' | sed "s/_1.$ext//");
   rv=$(echo $fw | sed "s/_1.$ext/_2.$ext/");
   #cutadapt -g "^GTGYCAGCMGCCGCGGTAA" --match-read-wildcards -O 19 -G "GGACTACNVGGGTWTCTAAT" -o cascabel_project/runs/Test/summer_data/demultiplexed/$sample.1_noprimer.fastq -p cascabel_project/runs/Test/summer_data/demultiplexed/$sample.2_noprimer.fastq $file $rv > cascabel_project/runs/Test/summer_data/demultiplexed/$sample.cutadapt.log
   cutadapt -g "$3" -G "$4" --match-read-wildcards -O $5 $6 -o $1/$sample.1_noprimer.fastq.gz -p $1/$sample.2_noprimer.fastq.gz $fw $rv > $1/$sample.cutadapt.log
   mv $1/$sample.1_noprimer.fastq.gz $fw
   mv $1/$sample.2_noprimer.fastq.gz $rv
   res=$(grep "(passing filters)" $1/$sample.cutadapt.log | awk '{print $5"\t"$6}')
   echo -e $sample"\t"$res >> $1/summary.txt
done

