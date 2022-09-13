"""
This Script use a mapping file to identify primer sequences and run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess
from benchmark_utils import countFasta
from benchmark_utils import countFastaGZ

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
with open(snakemake.input[0]) as mappingFile:
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
                if col == "ReverseLinkerPrimerSequence"  or col == "RvLinkerPrimerSequence" or col == "ReversePrimer"  :
                    idx_rv_primer=c
                elif col == "LinkerPrimerSequence":
                    idx_fw_primer=c
                elif col == "ReverseLinkerPrimerSequenceRevCom"  or col == "RvLinkerPrimerSequenceRevCom" or col == "ReversePrimerRevCom":
                    idx_rv_revcomp_primer=c                   
                c=c+1
            #if there is no "ReverseLinkerPrimerSequence" we look for the ReverseLinkerPrimerSequenceRevCom
            if idx_rv_primer == -1 and idx_rv_revcomp_primer !=1:
                idx_rv_primer=idx_rv_revcomp_primer
                isRC = True 
        elif not line.startswith("#"):
            if idx_rv_primer != -1:
                #if the valuee is the reverse complemented now we want the 5' to 3' orientation so rev com again.
                if isRC:
                    primer_by_sample[columns[0]]=[columns[idx_fw_primer],reverse_complement(columns[idx_rv_primer])]
                else:
                    primer_by_sample[columns[0]]=[columns[idx_fw_primer],columns[idx_rv_primer]]
            elif idx_fw_primer != -1:
                primer_by_sample[columns[0]]=[columns[idx_fw_primer]]
            else:
                print("\033[91m ERROR: LinkerPrimerSequence not found on mapping file: "+ snakemake.input[0] +" \033[0m")
                exit(1)


# List files
fq_files = [f for f in os.listdir(snakemake.params[0]) if f.endswith("_1."+snakemake.params[2])]
if not os.path.exists(snakemake.params[0]+"/reads_discarded_primer/"):
    os.makedirs(snakemake.params[0]+"/reads_discarded_primer/")
if not os.path.exists(snakemake.params[0]+"/primer_removed/"):
    os.makedirs(snakemake.params[0]+"/primer_removed/")
if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files")
summ_file = open(snakemake.output[0],"w")
summ_file2 = open(snakemake.params[4],"w")
summ_file.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
summ_file2.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
for fw in fq_files:
    sample=fw.replace("_1."+snakemake.params[2],"")
    fw_fq= snakemake.params[0]+"/"+fw
    rv=fw.replace("_1."+snakemake.params[2],"_2."+snakemake.params[2])
    rv_fq= snakemake.params[0]+"/"+rv
    if sample in primer_by_sample:
        runCutAdapt=False
        discard_untrimmed=""
        extra_params=snakemake.params[1] 
        if len(primer_by_sample[sample])>1 and snakemake.params[3] == "PE":
            if "--discard-untrimmed" in snakemake.params[1]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz --untrimmed-paired-output  "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_2.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","")
            #print("cutadapt -g "+ primer_by_sample[sample][0] + " -G " + primer_by_sample[sample][1] + " " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log")
            subprocess.run(["cutadapt -g "+ primer_by_sample[sample][0] + " -G " + primer_by_sample[sample][1] + " " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)
            runCutAdapt=True
            #subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
            #subprocess.run( ["cutadapt "+ primer_set +" "+snakemake.params[0]+" -o "+snakemake.output[0] + " " + snakemake.input[0]+ ">"+ snakemake.output[1]],stdout=subprocess.PIPE, shell=True)
        elif len(primer_by_sample[sample])>=1 and snakemake.params[3] == "SE":
            if "--discard-untrimmed" in snakemake.params[0]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","") 
            subprocess.run(["cutadapt -g "+ primer_by_sample[sample][0] +" " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz "+ discard_untrimmed + " " + fw_fq + " >> "+ snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)  
            #subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
            runCutAdapt=True
        elif len(primer_by_sample[sample])==1 and snakemake.params[3] == "PE":
            print("\033[91m ERROR: Found forward and reverse reads, but only one primer was supplied \033[0m")
            print("sample: "+sample + " primer " + primer_by_sample[sample][0])
            summ_file.close()
            summ_file2.close()
            exit(1)
        
        if runCutAdapt:
            if snakemake.params[2].endswith("gz"):
                reads_ori=countFastaGZ(fw_fq,True)
                reads_after=countFastaGZ(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz",True)
            else:
                reads_ori=countFasta(fw_fq,True)
                reads_after=countFasta(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq",True)
            prcOK="{:.2f}".format(float((reads_after/reads_ori)*100))
            summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
            summ_file2.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
 
    else:
        print("\033[93m WARNING: No primers found for sample: "+sample +" \033[0m")
        summ_file.close()
        summ_file2.close()
        exit(1)
summ_file.close() 
summ_file2.close()

