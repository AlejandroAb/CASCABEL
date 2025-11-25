#!/bin/bash
#This bash script concatenates the run of 2 cutadapt 
#Params:
#$1 directory with fastq files
#$2 -t Specify  the  minimum length of reads after trimming the low quality part
#$3 -v The minimum overlap size
#$4 -p p-value for the statistical test
#$5 -j Number of threads to use
#$6 -q  quality  score  threshold for trimming the low quality  part  of  a  read
#$7 output path
#$8 only some times, extra_params


for file in $1/*_1.fastq.gz ; 
  do
   id=$(echo $file | awk -F"/" '{print $NF}' | sed 's/_1.fastq.gz//')
   rv=$1/${id}_2.fastq.gz;
   echo "Processing: " $id 
   if  [ $# -lt 8 ]; then
     echo " pear -f "$file " -r "${rv}" -t "$2" -v "$3" -p "$4" -j "$5" -q "$6" -o " $7"/"${id} ">" $7"/"${id}".log"
     pear -f $file -r ${rv} -t $2 -v $3 -p $4 -j $5 -q $6 -o $7/${id} > $7/${id}.log; 
   else
     echo " pear -f "$file " -r "${rv}" -t "$2" -v "$3" -p "$4" -j "$5" -q "$6 "  " $8" -o " $7"/"${id} ">" $7"/"${id}".log"
     pear -f $file -r ${rv} -t $2 -v $3 -p $4 -j $5 -q $6 $8 -o $7/${id} > $7/${id}.log;
   fi 
   echo "done"
done
