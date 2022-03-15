#!/bin/bash
#need 6 arguments
#$1 path to demultiplexed files
#$2 extension fastq.gz or only fastq
#$3 FW primer (for read $sample_1.fastq.gz)
#$4 RV primer (for read $sample_2.fastq.gz)
#$5 Require MINLENGTH overlap between read and adapter for an adapter to be found
#$6 extra params for cutadapt
mkdir -p $1/primer_removed 
for fw in $1/*_1.$2 ; do
   ext=$2
   sample=$(echo $fw | awk -F"/" '{print $NF}' | sed "s/_1.$ext//");
   rv=$(echo $fw | sed "s/_1.$ext/_2.$ext/");
   cutadapt -g "$3" -G "$4" --match-read-wildcards -O $5 $6 -o $1/primer_removed/${sample}_1.fastq.gz -p $1/primer_removed/${sample}_2.fastq.gz $fw $rv > $1/primer_removed/$sample.cutadapt.log
   res=$(grep "(passing filters)" $1/primer_removed/$sample.cutadapt.log | awk '{print $5"\t"$6}')
   echo -e $sample"\t"$res >> $1/summary.txt
done

