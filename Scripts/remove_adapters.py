"""
This Script use a mapping file to identify primer sequences and run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess
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
if snakemake.config["cutAdapters"].lower() == "metadata":
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
                    if col == "ReversePrimer" or col == "LinkerPrimerSequenceReverse"  or col == "ReverseLinkerPrimerSequence"  or col == "RvLinkerPrimerSequence" :
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
                    #here, we were creating a dic with sample:primer, but is not needed anymore  
                    #primer_by_sample[columns[0]]=[columns[idx_fw_primer],columns[idx_rv_primer]]
                    #for primer in uniq_primers:
                    if columns[idx_fw_primer]+columns[idx_rv_primer] not in uniq_primers:
                        if isRC:
                            uniq_primers[columns[idx_fw_primer]+columns[idx_rv_primer]]=[columns[idx_fw_primer],columns[idx_rv_primer]]
                        else:
                            uniq_primers[columns[idx_fw_primer]+columns[idx_rv_primer]]=[columns[idx_fw_primer],reverse_complement(columns[idx_rv_primer])]    
                else:
                    #primer_by_sample[columns[0]]=[columns[idx_fw_primer]]
                    if columns[idx_fw_primer] not in uniq_primers:
                        uniq_primers[columns[idx_fw_primer]]=[columns[idx_fw_primer]]


#run_by_primer = True
i=0
primer_set = ""
no_primer=""
extra=snakemake.params[0]

if "--discard-untrimmed" in snakemake.params[0]:
    no_primer=" --untrimmed-output " + snakemake.params[2]
    extra=snakemake.params[0].replace("--discard-untrimmed","")

if snakemake.config["cutAdapters"].lower() == "metadata":
#if run_by_primer :
    for key in uniq_primers:
        i=i+1
        if len(uniq_primers[key])>1:
            #Now we always store the revcom of the reverse primer 
            #rc=reverse_complement(uniq_primers[key][1])
            #Here we use to build all the combinations, but maybe we only want the FWP...rc(RVP)
            #primer_set=primer_set+" -a "+uniq_primers[key][0]+"..."+uniq_primers[key][1]+" -a "+uniq_primers[key][0]+"..."+ rc + " " 
            primer_set=primer_set+" -g "+uniq_primers[key][0]+"..."+uniq_primers[key][1]+" "
            #print("cutadapt -f fasta -a "+uniq_primers[key][0]+"..."+uniq_primers[key][1]+" -a "+uniq_primers[key][0]+"..."+ rc + " "+snakemake.params[0]+" -o "+snakemake.output[0]+"_"+str(i) +" "+snakemake.input[0] + " > "+ snakemake.output[1]+"_"+str(i))
            #subprocess.run( ["cutadapt -f fasta -a "+uniq_primers[key][0]+"..."+uniq_primers[key][1]+" -a "+uniq_primers[key][0]+"..."+ rc + " "+snakemake.params[0]+" -o "+snakemake.output[0]+"_"+str(i) + " " + snakemake.input[0]+ ">"+ snakemake.output[1]+"_"+str(i)],stdout=subprocess.PIPE, shell=True)
        else:
            #rc=reverse_complement(uniq_primers[key][0])
            #primer_set=primer_set+" -a "+uniq_primers[key][0]+" -a " + rc
            primer_set=primer_set+" -g "+uniq_primers[key][0] +" "
    with open(snakemake.params[1], "w") as primers:
            primers.write(primer_set)
            primers.close()

    subprocess.run( ["cutadapt -f fasta "+ primer_set +" "+extra+" -o "+snakemake.output[0] + " "+ no_primer +" " + snakemake.input[0]+ ">"+ snakemake.output[1]],stdout=subprocess.PIPE, shell=True)
else:
    subprocess.run( ["cutadapt -f fasta "+ snakemake.config["cutadapt"]["adapters"] +" "+extra+" -o "+snakemake.output[0] + " "+ no_primer +" " + snakemake.input[0]+ ">"+ snakemake.output[1]],stdout=subprocess.PIPE, shell=True)

if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files")
subprocess.run( ["cat "+ snakemake.output[0]+"| grep '^>' | cut -f1 -d' ' | sed 's/>// ; s/_[0-9]*$//' |  sort | uniq -c | awk '{print $2\"\\t\"$1}' > " + snakemake.params[3]+".tmp1"],stdout=subprocess.PIPE, shell=True)
subprocess.run( ["cat "+ snakemake.input[0]+"| grep '^>' | cut -f1 -d' ' | sed 's/>// ; s/_[0-9]*$//' |  sort | uniq -c | awk '{print $2\"\\t\"$1}'| awk -F'\t' 'NR==FNR{h[$1]=$2;next} BEGIN{print \"Sample\\tReads_before_cutadapt\\tSurviving_reads\\tPrc_surviving_reads\"}{if(h[$1]){print $1\"\\t\"h[$1]\"\\t\"$2\"\\t\"($2/h[$1])*100\"%\"}else{print $1\"\\t\"$2\"\\t0\\t0%\"}}' - "+snakemake.params[3]+".tmp1 > "+ snakemake.params[3]],stdout=subprocess.PIPE, shell=True)
os.remove(snakemake.params[3]+".tmp1")
#sample=snakemake.params[3].split("/")[3].split("_")[0]
#summ_file = open(snakemake.params[3],"w")
#summ_file.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
#reads_ori=countFasta(snakemake.input[0],False)
#reads_after=countFasta(snakemake.output[0],False)
#prcOK="{:.2f}".format(float((reads_after/reads_ori)*100))
#summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
#summ_file.close()    


