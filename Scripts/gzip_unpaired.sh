#!/bin/bash
#This script compress the unpaired fastq files fom pear and place them in
#the demultiplexed/unpaired folder for a potential ASV unpaired workflow.
#Params:
#$1 input directory
#$2 output directory
#$3 output log/summary file

echo -e "#SampeID\tAssembled_reads\tDiscarded_reads\tNot_assembled_reads" > $3

for file in $1/*.unassembled.forward.fastq; 
do
 id=$(echo $file | awk -F"/" '{print $NF}' | sed 's/.unassembled.forward.fastq//');
 rv=$(echo $file | awk -F"/" '{print $NF}' | sed 's/.unassembled.forward.fastq/.unassembled.reverse.fastq/');
 cat $1/${id}.log |  awk -F":" -v id=${id} '$0 ~ "^Assembled reads"{if(!assembled)assembled = $2} 
  		                                    $0 ~ "^Discarded reads"{if(!discard)discard = $2} 
  		                                    $0 ~ "^Not assembled reads"{not_ass = $2} 
  		                                    END{printf "%s\t%s\t%s\t%s\n",id,assembled,discard,not_ass}' >> $3
 
 cat $file |  gzip -9 > $2/${id}_1.fastq.gz
 cat $1/$rv |  gzip -9 > $2/${id}_2.fastq.gz

 
done

 
