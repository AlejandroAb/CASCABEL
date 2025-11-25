#!/bin/bash
#This bash script concatenates the run of 2 cutadapt 
#Params:
#$1 directory with fastq files
#$2 output directory for concatenated fastq files


for file in $1/*_1.fastq.gz ; 
  do
   id=$(echo $file | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//')
   rv_fbc=$1/${id}_fw.fastq.gz
   rv_rbc=$(echo $file | sed 's/_1.fastq.gz/_2.fastq.gz/')
   if [ -f $rv_fbc ]; then
     #with legend on the 3rd row
     #cat $file <(zcat ${rv_fbc} | awk '{if(NR%4==1){gsub(" 2:N:"," 1:N:",$0)}else if(NR%4==3){gsub("+","+ from_read_2",$0)}; print $0}' | gzip -9 ) > $2/${id}_1.fastq.gz
     cat $file <(zcat ${rv_fbc} | awk '{if(NR%4==1){gsub(" 2:N:"," 1:N:",$0)}; print $0}' | gzip -9 ) > $2/${id}_1.fastq.gz
     #with legend on the 3rd row
     #cat ${rv_rbc} <(zcat $1/${id}_rv.fastq.gz | awk '{if(NR%4==1){gsub(" 1:N:"," 2:N:",$0)}else if(NR%4==3){gsub("+","+ from_read_1",$0)}; print $0}' | gzip -9) > $2/${id}_2.fastq.gz
     cat ${rv_rbc} <(zcat $1/${id}_rv.fastq.gz | awk '{if(NR%4==1){gsub(" 1:N:"," 2:N:",$0)}; print $0}' | gzip -9) > $2/${id}_2.fastq.gz
     echo $id "FW" $file ${rv_fbc} 
     echo $id "RV" ${rv_rbc} $1/${id}_rv.fastq.gz
   else
     cp $file > $2/${id}_1.fastq.gz
     cp ${rv_rbc} > $2/${id}_2.fastq.gz
     echo $id "FW" $file "no_second_file"
     echo $id "RV" ${rv_rbc} "no_second_file" 
   fi
done
