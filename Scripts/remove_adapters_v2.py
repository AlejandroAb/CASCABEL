"""
This Script use a mapping file to identify primer sequences and run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess
from benchmark_utils import countFasta
from sys import stdin

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y','S':'S','W':'W','K':'M','M':'K','N':'N','B':'V','V':'B','D':'H','H':'D'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

primer_by_sample={}
uniq_primers={}
idx_fw_primer=-1   # default for qiime (col 3)
idx_rv_primer=-1   # new field 
idx_rv_revcomp_primer=-1
isRC = False  
primer_set = ""
no_primer = ""
extra=snakemake.params[0]
log_by_sample="Sample\tInitial reads\tSurviving reads\n"
if "--discard-untrimmed" in snakemake.params[0]:
    no_primer=" --untrimmed-output " + snakemake.params[2]
    extra=snakemake.params[0].replace("--discard-untrimmed","")

if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files")

if snakemake.config["primers"]["remove"].lower() == "metadata":
    with open(snakemake.input[1]) as mappingFile:
        l=0
        for line in mappingFile:
            l=l+1;
            columns = line.split('\t')
            #the header is always at row 1 and must contain these first 3 fields (qiime specs):
            #SampleID BarcodeSequence LinkerPrimerSequence Description
            if l==1 :
                c=0
                #Find target headers
                for col in columns:
                    if col == "ReversePrimer" or col == "LinkerPrimerSequenceReverse"  or col == "ReverseLinkerPrimerSequence"  or col == "RvLinkerPrimerSequence" or col == "ReversePrimerSequence" :
                        idx_rv_primer=c
                    elif col == "LinkerPrimerSequence":
                        idx_fw_primer=c
                    elif col == "ReverseLinkerPrimerSequenceRevCom"  or col == "ReversePrimerRevCom":
                        idx_rv_revcomp_primer=c
                        isRC=True
                    c=c+1
                if isRC:
                    idx_rv_primer=idx_rv_revcomp_primer 
            elif not line.startswith("#"):
                if idx_rv_primer != -1:
                    #here, we are creating a dic with sample:primer
                    if isRC:  
                        primer_by_sample[columns[0]]=[columns[idx_fw_primer],columns[idx_rv_primer]]
                    else:
                        primer_by_sample[columns[0]]=[columns[idx_fw_primer],reverse_complement(columns[idx_rv_primer])]
                    #for primer in uniq_primers:
                    if columns[idx_fw_primer]+columns[idx_rv_primer] not in uniq_primers:
                        if isRC:
                            uniq_primers[columns[idx_fw_primer]+columns[idx_rv_primer]]=[columns[idx_fw_primer],columns[idx_rv_primer]]
                        else:
                            uniq_primers[columns[idx_fw_primer]+columns[idx_rv_primer]]=[columns[idx_fw_primer],reverse_complement(columns[idx_rv_primer])]    
                else:
                    primer_by_sample[columns[0]]=[columns[idx_fw_primer]]
                    if columns[idx_fw_primer] not in uniq_primers:
                        uniq_primers[columns[idx_fw_primer]]=[columns[idx_fw_primer]]
        mappingFile.close()
                         
    #If we have more than one different pair of primers, we run cutadapt by sample
    #otherwise we run only one instance
    if len(uniq_primers) >1:
        #create tmp dir
        if not os.path.exists(snakemake.params[4]):
            os.makedirs(snakemake.params[4])
        else: #it exists and most lickly we want to delete all its content.
            subprocess.run( ["rm -fr " + snakemake.params[4]+"*"],stdout=subprocess.PIPE, shell=True)
        #split the reads
        #If we are running this, it comes from our demultiplexing, and thus we have fasta headers like this:
        #><sample>_###  so we remove the _###
        subprocess.run(["cat "+ snakemake.input[0]+ " |  awk '{if($0 ~ \"^>\"){sample=$1; header=$0; gsub(\">\",\"\",sample);gsub(\"_[0-9].*\",\"\",sample);}else{print header\"\\n\"$0 >> \""+snakemake.params[4]+"\"sample\".fasta\"} }'"],stdout=subprocess.PIPE, shell=True)
        all_primers=""
        for file in os.listdir(snakemake.params[4]):
            #file only has the name of the file, the path is already discarded
            #the function os.path.splitext strip the extension
            sample=os.path.splitext(file)[0]
            no_primer=""
            extra=""
            if "--discard-untrimmed" in snakemake.params[0]:
                no_primer=" --untrimmed-output " + snakemake.params[4]+sample+"_untrimmed.fasta"
                extra=snakemake.params[0].replace("--discard-untrimmed","")
            tmp_out = snakemake.params[4]+sample+"_trimmed.fasta"
            tmp_log = snakemake.params[4]+sample+".log"
            if sample in primer_by_sample:
                if len(primer_by_sample[sample])>1:
                    primer_set=" -g "+primer_by_sample[sample][0]+"..."+primer_by_sample[sample][1]+" "
                else:
                    primer_set=" -g "+primer_by_sample[sample][0]
                #run cutadapt by sample
                subprocess.run(["echo \"Processing sample\" " + sample + "\n >> "+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
                subprocess.run( ["cutadapt  "+ primer_set +" "+extra+" -o "+tmp_out + " "+ no_primer +" " + snakemake.params[4]+file+ ">>"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
                #stats by sample
                initialReads=countFasta(snakemake.params[4]+file,False)
                survivingReads=countFasta(tmp_out,False)
                prc = "{:.2f}".format(float((survivingReads/initialReads)*100))
                log_by_sample=log_by_sample+sample+"\t"+str(initialReads)+"\t"+str(survivingReads)+" ("+prc+"%)\n"
                all_primers=all_primers+sample+"\t"+primer_set+"\n"
            else:
                print("\033[91mNo primers found for sample:"+ sample+ " \033[0m")
                print("\033[91mPlease make sure to have the sample included in the mapping file: "+snakemake.input[1]+"  \033[0m")
                print("\033[91mAborting the pipeline \033[0m")
                exit(1)

        #merge results
        subprocess.run( ["cat  "+ snakemake.params[4]+"*_trimmed.fasta > "+ snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
        with open(snakemake.params[1], "a") as primers:
            primers.write(all_primers)
            primers.close()
        #subprocess.run( ["cat  "+ snakemake.params[4]+"*_untrimmed.fasta > " snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
    else: #only run one cutadapt instance
        new_key = list(uniq_primers)
        if len(uniq_primers[new_key[0]])>1: #is PE?
            primer_set=" -g "+uniq_primers[new_key[0]][0]+"..."+uniq_primers[new_key[0]][1]+" "
        else: #is SE
            primer_set=" -g "+uniq_primers[new_key[0]][0]
        subprocess.run( ["cutadapt  "+ primer_set +" "+extra+" -o "+snakemake.output[0] + " "+ no_primer +" " + snakemake.input[0]+ ">"+  snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
        with open(snakemake.params[1], "a") as primers:
            primers.write(primer_set)
            primers.close()
else: #values come at the CFG, run only once
    primer_set="-g " + snakemake.config["primers"]["fw_primer"]
    if len(snakemake.config["primers"]["rv_primer"]) > 2 :
        primer_set=primer_set+"..."+reverse_complement(snakemake.config["primers"]["rv_primer"])

    subprocess.run( ["cutadapt  "+ primer_set  +" "+extra+" -o "+snakemake.output[0] + " "+ no_primer +" " + snakemake.input[0]+ ">"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
  #  primer_set = snakemake.config["cutadapt"]["adapters"]
    with open(snakemake.params[1], "w") as primers:
        primers.write(primer_set)
        primers.close()

initialReads=countFasta(snakemake.input[0],False)
survivingReads=countFasta(snakemake.output[0],False)
prc=float((survivingReads/initialReads)*100)
prc_str = "{:.2f}".format(float((survivingReads/initialReads)*100))

user_input="0"
while (user_input != "1" and user_input !=  "2"):
    print("\033[91m This step removes primers \033[0m")
    print("\033[93m Total number of initial reads: " + str(initialReads) + " \033[0m")
    print("\033[93m Total number of surviving reads: " + str(survivingReads) + " = "+ prc_str + "% \033[0m")
    print("\033[93m You can find cutadapt's log file at: " + snakemake.params[5] +"\n \033[0m")
    if snakemake.config["interactive"] != "F" or prc < snakemake.config["primers"]["min_prc"]:
        print("\033[92m What would you like to do? \033[0m")
        print("\033[92m 1. Continue with the workflow. \033[0m")
        print("\033[92m 2. Interrupt the workflow. \033[0m")
        if snakemake.config["primers"]["remove"].lower() == "metadata" and  len(uniq_primers)>1:
            print("\033[92m 3. Print results by sample. \033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = user_input[:-1]
        if user_input == "2":
            print("\033[91m Aborting workflow... \033[0m")
            #delete target outpu (snakemake also does it)
            subprocess.run( ["rm -f "+ snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
            #delete cutadapt mp directory
            subprocess.run( ["rm -fr " + snakemake.params[4]],stdout=subprocess.PIPE, shell=True)
            #delete all the concatenated log files
            subprocess.run( ["rm -f " + snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
            #delete primers file
            subprocess.run( ["rm -f " + snakemake.params[1]],stdout=subprocess.PIPE, shell=True)
            exit(1)
        if user_input == "3":
            print(log_by_sample)
    else:
        print("\033[93m" +" Interactive mode off \033[0m")
        print("\033[93m" +" Removing primers...\033[0m")
        user_input="1"
# if we ran multiple cutadap tasks, now delete tmp files and logs. 
if snakemake.config["primers"]["remove"].lower() == "metadata" and  len(uniq_primers)>1: 
    print("\033[96mCleaning intermediate files...\033[0m")
    subprocess.run( ["rm -fr " + snakemake.params[4]],stdout=subprocess.PIPE, shell=True)

#Summarize results    
subprocess.run( ["cat "+ snakemake.output[0]+"| grep '^>' | cut -f1 -d' ' | sed 's/>// ; s/_[0-9]*$//' |  sort | uniq -c | awk '{print $2\"\\t\"$1}' > " + snakemake.params[3]+".tmp1"],stdout=subprocess.PIPE, shell=True)
subprocess.run( ["cat "+ snakemake.input[0]+"| grep '^>' | cut -f1 -d' ' | sed 's/>// ; s/_[0-9]*$//' |  sort | uniq -c | awk '{print $2\"\\t\"$1}'| awk -F'\t' 'NR==FNR{h[$1]=$2;next} BEGIN{print \"Sample\\tReads_before_cutadapt\\tSurviving_reads\\tPrc_surviving_reads\"}{if(h[$1]){print $1\"\\t\"h[$1]\"\\t\"$2\"\\t\"($2/h[$1])*100\"%\"}else{print $1\"\\t\"$2\"\\t0\\t0%\"}}' - "+snakemake.params[3]+".tmp1 > "+ snakemake.params[3]],stdout=subprocess.PIPE, shell=True)
os.remove(snakemake.params[3]+".tmp1")
exit(0)
#sample=snakemake.params[3].split("/")[3].split("_")[0]
#summ_file = open(snakemake.params[3],"w")
#summ_file.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
#reads_ori=countFasta(snakemake.input[0],False)
#reads_after=countFasta(snakemake.output[0],False)
#prcOK="{:.2f}".format(float((reads_after/reads_ori)*100))
#summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
#summ_file.close()    


