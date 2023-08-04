"""
This Script use the primer values from the config file to run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess
from sys import stdin
from benchmark_utils import countFasta
from benchmark_utils import countFastaGZ
import shutil

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y','S':'S','W':'W','K':'M','M':'K','N':'N','B':'V','V':'B','D':'H','H':'D'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

# List files
fq_files = [f for f in os.listdir(snakemake.params[0]) if f.endswith("_1."+snakemake.params[2])]
if not os.path.exists(snakemake.params[0]+"/reads_discarded_primer/") and "--discard-untrimmed" in snakemake.params[1]:
    os.makedirs(snakemake.params[0]+"/reads_discarded_primer/")
if not os.path.exists(snakemake.params[0]+"/primer_removed/"):
    os.makedirs(snakemake.params[0]+"/primer_removed/")
if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files")

summ_file = open(snakemake.output[0],"w") # this iss a log for the wf
summ_file2 = open(snakemake.params[4],"w") # this is for the report
summ_file.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n");
summ_file2.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n");
log_str = "Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n"
log_zero = "Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n"
has_zero_length_reads = False
zero_samples = 0;
to_remove = []

for fw in fq_files:
    sample=fw.replace("_1."+snakemake.params[2],"")
    fw_fq= snakemake.params[0]+"/"+fw
    rv=fw.replace("_1."+snakemake.params[2],"_2."+snakemake.params[2])
    rv_fq= snakemake.params[0]+"/"+rv
    discard_untrimmed=""
    extra_params=snakemake.params[1]
#Count reads before trimming
    if snakemake.params[2].endswith("gz"):
        reads_ori=countFastaGZ(fw_fq,True)
    else:
        reads_ori=countFasta(fw_fq,True)
#no cutadapt if no reads 
    if reads_ori > 0:
        if snakemake.params[3] == "PE":
            if "--discard-untrimmed" in snakemake.params[1]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz --untrimmed-paired-output  "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_2.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","")
            subprocess.run(["cutadapt -g "+ snakemake.config["primers"]["fw_primer"]  + " -G " + snakemake.config["primers"]["rv_primer"]  + " " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"]+" -m "+ snakemake.config["primers"]["min_length"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)
        elif snakemake.params[3] == "SE":
            if "--discard-untrimmed" in snakemake.params[0]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","") 
            subprocess.run(["cutadapt -g "+ snakemake.config["primers"]["fw_primer"] +" " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"]+" -m "+ snakemake.config["primers"]["min_length"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz "+ discard_untrimmed + " " + fw_fq + " >> "+ snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)  
        
        if snakemake.params[2].endswith("gz"):
            reads_after=countFastaGZ(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz",True)
        else:
            reads_after=countFasta(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq",True)
        
        prcOK="{:.2f}".format(float((reads_after/reads_ori)*100))

    else:
        reads_after = 0
        prcOK="{:.2f}".format(float((reads_after/1)*100))
        to_copy=snakemake.params[0]+"/primer_removed/"+sample+"_1."+snakemake.params[2]
        os.symlink(fw_fq,to_copy)
        if snakemake.params[3] == "PE":
            to_copy_rv=snakemake.params[0]+"/primer_removed/"+sample+"_2."+snakemake.params[2]
            os.symlink(rv_fq,to_copy_rv)

    if reads_after < 1:
        has_zero_length_reads = True
        log_zero = log_zero + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
        to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_1."+snakemake.params[2])
        if snakemake.params[3] == "PE":
            to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_2."+snakemake.params[2])
        zero_samples = zero_samples + 1

    log_str = log_str + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
    summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
    summ_file2.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");

summ_file.close()
summ_file2.close()

user_input="0"
show_menu = True
if zero_samples > 0:
    while show_menu:
        print("\033[91m\n###########  Primer removal validation    ###########\033[0m")
        print("\033[91m You have " + str(zero_samples) + " samples without reads surviving filters. \033[0m")
        print("\033[92m LIBRARY: "+snakemake.wildcards.sample+" \033[0m")
        print("\033[92m cutadapt_log: "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log \033[0m")
        print("\033[93m Please select one of the following options: \033[0m")
        print("\033[93m   1. Print samples with 0 reads \033[0m")
        print("\033[93m   2. Print summary (all the samples) \033[0m")
        print("\033[93m   3. Remove from this analysis samples with 0 reads\033[0m")
        print("\033[93m      and continue with the workflow. \033[0m")
        print("\033[93m   4. Interrupt the workflow and re-do primer removal step. \033[0m")
        print("\033[93m      Adjust primer values in your configuration and/or mapping file \033[0m")
        print("\033[93m      and restart the pipeline. \033[0m")
        print("\033[93m      This action will remove:"+snakemake.params[0]+"/primer_removed \033[0m")
        print("\033[93m   5. Interrupt the workflow \033[0m")
        print("\033[93m   6. Continue with the workflow\n      (an error will be raised during dada2)\n      Pointless option... \033[0m")
        print("\033[93m Select an option: \033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = user_input[:-1]
        if user_input == "1":
            print(log_zero)
        elif user_input == "2":
            print(log_str)
        elif user_input == "3":
            for file in to_remove:
                newn = file+"_NOK"
                os.rename(file, newn)
                show_menu = False
        elif user_input == "4":
            shutil.rmtree(snakemake.params[0]+"/primer_removed")
            exit(1) 
        elif user_input == "5":
            exit(1)

exit(0)
