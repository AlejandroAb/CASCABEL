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
from sys import stdin
import shutil

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y','S':'S','W':'W','K':'M','M':'K','N':'N','B':'V','V':'B','D':'H','H':'D'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

#from Bio.Seq import Seq

run=snakemake.params[5]
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
                if col == "ReverseLinkerPrimerSequence"  or col == "RvLinkerPrimerSequence" or col == "ReversePrimer" or col == "ReversePrimerSequence"  :
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
if not os.path.exists(snakemake.wildcards.PROJECT+"/runs/"+run+"/report_files"):
    os.makedirs(snakemake.wildcards.PROJECT+"/runs/"+run+"/report_files")
summ_file = open(snakemake.output[0],"w")
summ_file2 = open(snakemake.params[4],"w")
summ_file.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
summ_file2.write("Sample\tReads_before_cutadapt\tSurviving_reads\tPrc_surviving_reads\n")
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
    #Count reads before trimming
    if snakemake.params[2].endswith("gz"):
        reads_ori=countFastaGZ(fw_fq,True)
    else:
        reads_ori=countFasta(fw_fq,True)
    #no cutadapt if no reads
    #if reads_ori > 0:

    if sample in primer_by_sample:
        runCutAdapt=False
        discard_untrimmed=""
        extra_params=snakemake.params[1] 
        if len(primer_by_sample[sample])>1 and reads_ori < 1:
            reads_after = 0
            prcOK="{:.2f}".format(float((reads_after/1)*100))
            to_copy=snakemake.params[0]+"/primer_removed/"+sample+"_1."+snakemake.params[2]
            os.symlink(fw_fq,to_copy)
            if snakemake.params[3] == "PE":
                to_copy_rv=snakemake.params[0]+"/primer_removed/"+sample+"_2."+snakemake.params[2]
                os.symlink(rv_fq,to_copy_rv)
            runCutAdapt=True
    
        elif len(primer_by_sample[sample])>1 and snakemake.params[3] == "PE" :
            if "--discard-untrimmed" in snakemake.params[1]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz --untrimmed-paired-output  "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_2.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","")
            #print("cutadapt -g "+ primer_by_sample[sample][0] + " -G " + primer_by_sample[sample][1] + " " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log")
            subprocess.run(["cutadapt -g "+ primer_by_sample[sample][0] + " -G " + primer_by_sample[sample][1]+" -m "+ snakemake.config["primers"]["min_length"]  + " " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)
            runCutAdapt=True
            #subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
            #subprocess.run( ["cutadapt "+ primer_set +" "+snakemake.params[0]+" -o "+snakemake.output[0] + " " + snakemake.input[0]+ ">"+ snakemake.output[1]],stdout=subprocess.PIPE, shell=True)
        elif len(primer_by_sample[sample])>=1 and snakemake.params[3] == "SE":
            if "--discard-untrimmed" in snakemake.params[0]:
                discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz"
                extra_params=snakemake.params[1].replace("--discard-untrimmed","") 
            subprocess.run(["cutadapt -g "+ primer_by_sample[sample][0] +" -m "+ snakemake.config["primers"]["min_length"] +" " +extra_params+" -O "+ snakemake.config["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz "+ discard_untrimmed + " " + fw_fq + " >> "+ snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)  
            #subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
            runCutAdapt=True
        elif len(primer_by_sample[sample])==1 and snakemake.params[3] == "PE":
            print("\033[91m ERROR: Found forward and reverse reads, but only one primer was supplied \033[0m")
            print("sample: "+sample + " primer " + primer_by_sample[sample][0])
            summ_file.close()
            summ_file2.close()
            exit(1)
        
        if runCutAdapt and reads_ori > 0:
            if snakemake.params[2].endswith("gz"):
                reads_ori=countFastaGZ(fw_fq,True)
                reads_after=countFastaGZ(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz",True)
            else:
                reads_ori=countFasta(fw_fq,True)
                reads_after=countFasta(snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq",True)
            prcOK="{:.2f}".format(float((reads_after/reads_ori)*100))
            summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
            summ_file2.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
            log_str = log_str + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
            if reads_after < 1:
                has_zero_length_reads = True
                log_zero = log_zero + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
                to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_1."+snakemake.params[2])
                if snakemake.params[3] == "PE":
                    to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_2."+snakemake.params[2])
                zero_samples = zero_samples + 1


        elif runCutAdapt and reads_ori < 1:
            summ_file.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
            summ_file2.write(sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n");
            log_str = log_str + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
            log_zero = log_zero + sample+"\t"+str(reads_ori)+"\t"+str(reads_after)+"\t"+prcOK+"\n"
            zero_samples = zero_samples + 1
            has_zero_length_reads = True
            to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_1."+snakemake.params[2])
            if snakemake.params[3] == "PE":
                to_remove.append(snakemake.params[0]+"/primer_removed/"+sample+"_2."+snakemake.params[2])

 
    else:
        print("\033[93m WARNING: No primers found for sample: "+sample +" \033[0m")
        summ_file.close()
        summ_file2.close()
        exit(1)
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


