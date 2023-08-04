"""
This Script use a mapping file to identify primer sequences and run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess
from sys import stdin
#import benchmark_utils
from benchmark_utils import countFasta

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y','S':'S','W':'W','K':'M','M':'K','N':'N','B':'V','V':'B','D':'H','H':'D'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

#from Bio.Seq import Seq


primer_by_sample={}
uniq_primers={}
idx_fw_primer=-1   # default for qiime (col 3)
idx_rv_primer=-1  # new field 
idx_rv_revcomp_primer=-1
isRC = False
foundSample=False  
primer=""
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
            elif line.startswith(snakemake.params[4]):
                foundSample=True
                if idx_rv_primer != -1:
                    if isRC:
                        #fw_primer=columns[idx_fw_primer]
                        #rv_primer=columns[idx_rv_primer]
                        primer="-g "+columns[idx_fw_primer]+"..."+columns[idx_rv_primer]
                    else:
                        #fw_primer=columns[idx_fw_primer]
                        #rv_primer=reverse_complement(columns[idx_rv_primer])
                        primer="-g "+columns[idx_fw_primer]+"..."+reverse_complement(columns[idx_rv_primer])
                else:
                    #fw_primer=columns[idx_fw_primer]
                    primer="-g "+columns[idx_fw_primer]


    if not foundSample:
        print("\033[91m" +"No primers found for sample:"+ snakemake.params[4]+ " \033[0m")
        print("\033[91mPlease make sure to have the sample included in the mapping file: "+snakemake.input[1]+"  \033[0m")
        print("\033[91m Aborting the pipeline \033[0m")
        exit(1)

elif snakemake.config["primers"]["remove"].lower() == "cfg":
    primer="-g " + snakemake.config["primers"]["fw_primer"]
    if snakemake.config["primers"]["rv_primer"].len() > 2 :
        primer=primer+"..."+reverse_complement(snakemake.config["primers"]["rv_primer"]) 

    
discard = True
if "--discard-untrimmed" in snakemake.params[0]:
    extra=snakemake.params[0].replace("--discard-untrimmed","")
else: 
    extra=snakemake.params[0]
    discard = False
   
#This file will contain the untrimmed reads for the first pass
no_primer=" --untrimmed-output " + snakemake.params[2]+".tmp"

if snakemake.config["primers"]["remove"].lower() == "metadata":
    subprocess.run( ["cutadapt  "+ primer +" "+ extra+" -o "+snakemake.output[0] + ".1 "+ no_primer +" " + snakemake.input[0]+ ">"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
else:
    subprocess.run( ["cutadapt  "+ primer +" "+extra+" -o "+snakemake.output[0] + ".1 "+ no_primer +" " + snakemake.input[0]+ ">"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
#    primer=snakemake.config["cutadapt"]["adapters"]
#comment above line because we just add the primer generation in the elif above....


initialReads=countFasta(snakemake.input[0],False)
disscardedReads=countFasta(snakemake.params[2]+".tmp",False)

#The "extra" var returns to the original values in the sense that if the user wants to disscard reads
# this option will be present on the final cutadapt command 
extra=snakemake.params[0]
#if we disscarded reads
if disscardedReads>0:
    #reverse complement disscardedReads
    subprocess.run( ["vsearch --fastx_revcomp "+ snakemake.params[2]+".tmp  --fastaout "+ snakemake.params[2]+".tmp2"],stdout=subprocess.PIPE, shell=True)
    if snakemake.config["primers"]["remove"].lower() == "metadata":
        if discard:
        #Run cutadapt on disscarded reads
            subprocess.run( ["cutadapt  "+ primer +" "+ extra+" -o "+snakemake.output[0] + ".2 " + snakemake.params[2]+".tmp2"+ ">>"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
        else:
            print("Running second cutadapt")
            print("cutadapt  "+ primer +" "+ extra+" -o "+snakemake.output[0] + ".2 --untrimmed-output "+ snakemake.output[0] + ".3 " + snakemake.params[2]+".tmp2"+ ">>"+ snakemake.params[5])
            subprocess.run( ["cutadapt  "+ primer +" "+ extra+" -o "+snakemake.output[0] + ".2 --untrimmed-output "+ snakemake.output[0] + ".3 " + snakemake.params[2]+".tmp2"+ ">>"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
            #reverse complement untrimmed disscardedReads
            subprocess.run( ["vsearch --fastx_revcomp "+ snakemake.output[0]+".3  --fastaout "+ snakemake.params[2]+".tmp3"],stdout=subprocess.PIPE, shell=True)
    else:
        if discard:
            subprocess.run( ["cutadapt "+ primer  +" "+extra+" -o "+snakemake.output[0] + ".2 " + snakemake.params[2]+".tmp2 >>"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
        else:
            subprocess.run( ["cutadapt "+ primer  +" "+extra+" -o "+snakemake.output[0] + ".2 --untrimmed-output "+ snakemake.output[0] + ".3 "  + snakemake.params[2]+".tmp2 >>"+ snakemake.params[5]],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["vsearch --fastx_revcomp "+ snakemake.output[0]+".3  --fastaout "+ snakemake.params[2]+".tmp3"],stdout=subprocess.PIPE, shell=True)
    
    if discard:        
        #Concatenate results
        subprocess.run( ["cat "+snakemake.output[0] + ".1 "+ snakemake.output[0] + ".2 > "+ snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
        #remove intermediate files: disscarded reads first round, disscarded reads RC, accepted reads first round, accepted reads second round
        subprocess.run( ["rm -f "+ snakemake.params[2]+".tmp "+ snakemake.params[2]+".tmp2 "+snakemake.output[0] + ".1 "+ snakemake.output[0] + ".2"],stdout=subprocess.PIPE, shell=True)
    else:
        #Concatenate results
        subprocess.run( ["cat "+snakemake.output[0] + ".1 "+ snakemake.output[0] + ".2 " + snakemake.params[2] + ".tmp3  > " + snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
        #remove intermediate files: disscarded reads first round, disscarded reads RC, accepted reads first round, accepted reads second round
        #subprocess.run( ["rm -f "+ snakemake.params[2]+".tmp "+ snakemake.params[2]+".tmp2 "+snakemake.output[0] + ".1 "+ snakemake.output[0] + ".2 "+ snakemake.params[2]+".tmp3"],stdout=subprocess.PIPE, shell=True)
else: #no reads to evaluate just rename file
    print("No untrimmed output!!!!")
    subprocess.run( ["mv "+snakemake.output[0] + ".1  > "+ snakemake.output[0]],stdout=subprocess.PIPE, shell=True)

survivingReads=countFasta(snakemake.output[0],False)
prc = float((survivingReads/initialReads)*100)
prc_str = "{:.2f}".format(float((survivingReads/initialReads)*100))

with open(snakemake.params[1], "w") as primers:
        primers.write(primer)
        primers.close()

print("\033[91m This step removes primers \033[0m")
print("\033[93m Total number of initial reads: " + str(initialReads) + " \033[0m")
print("\033[93m Total number of surviving reads: " + str(survivingReads) + " = "+ prc_str + "% \033[0m")
print("\033[93m You can find cutadapt's log file at: " + snakemake.params[5] +"\n \033[0m")
if snakemake.config["interactive"] != "F" or prc < snakemake.config["primers"]["min_prc"]:
    print("\033[92m Do you want to continue?(y/n): \033[0m")
    user_input = stdin.readline() #READS A LINE
    user_input = user_input[:-1]
    if user_input.upper() == "N" or user_input.upper() == "NO":
        subprocess.run( ["rm -f "+ snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
        exit(1)
else:
    print("\033[93m" +" Interactive mode off \033[0m")
    print("\033[93m" +" Removing primers...\033[0m")


if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files")
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


