#!/bin/bash
#This script prepares the fastq files for the OTU workflow.
#creates a single fasta file and rename the headers
#Params:
#$1 input directory
#$2 output log file
#$3 fasta output sequences assembled
#$4 fasta output sequences unaseembled
#$5 char to pair (could be empty)

echo -e "#SampeID\tAssembled_reads\tDiscarded_reads\tNot_assembled_reads" > $2

for file in $1/*.assembled.fastq; 
do
 id=$(echo $file | awk -F"/" '{print $NF}' | sed 's/.assembled.fastq//');
 cat $1/${id}.log |  awk -F":" -v id=${id} '$0 ~ "^Assembled reads"{if(!assembled)assembled = $2} 
  		                                    $0 ~ "^Discarded reads"{if(!discard)discard = $2} 
  		                                    $0 ~ "^Not assembled reads"{not_ass = $2} 
  		                                    END{printf "%s\t%s\t%s\t%s\n",id,assembled,discard,not_ass}' >> $2
 
 cat $file | awk -v id=${id} -v num=0 'NR%4==1{num++;print ">"id"_"num" "$0;} 
                                       NR%4==2{print $0}' >> $3;

 

 if  [ $# -lt 5 ]; then
    cat $1/${id}.unassembled.forward.fastq  | awk -v id=${id}  -v num=0 ' NR==FNR{if(NR%4==1){h=$1}else if(NR%4==2){seq[h]=$0}next;} {if(NR%4==1){h=$1;num++;}else if(NR%4==2){print ">"id"_un_"num" "h"\n"seq[h]$0}}' - $1/${id}.unassembled.reverse.fastq >>  $4
else
	cat $1/${id}.unassembled.forward.fastq  | awk -v id=${id}  -v num=0 -v ch="$5" 'NR==FNR{if(NR%4==1){h=$1}else if(NR%4==2){seq[h]=$0}next;}{if(NR%4==1){h=$1;num++}else if(NR%4==2){print ">"id"_un_"num"\n"seq[h] ch $0}}' - $1/${id}.unassembled.reverse.fastq >>  $4
    #this was for fastq using quality char to pair, for the moment we only use fasta
	#cat $1/${id}.unassembled.forward.fastq  | awk -v ch="$5" 'BEGIN{for(i=1;i<=length(ch);i++){chq=chq "G"}} NR==FNR{if(NR%4==1){h=$1}else if(NR%4==2){seq[h]=$0}next;}{if(NR%4==1){h=$1}else if(NR%4==2){print ">"h"\n"seq[h] ch $0}}' - $1/${id}.unassembled.reverse.fastq > 

fi
done

 
