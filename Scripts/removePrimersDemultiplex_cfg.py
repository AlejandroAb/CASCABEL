"""
This Script use a mapping file to identify primer sequences and run cutadapt.
If only found "LinkerPrimerSequence" runs cutadapt with option -g LPS
If found "ReverseLinkerPrimerSequence" runs cutadapt with option -a LPS...rc(RvLPS)
If found "ReverseLinkerPrimerSequenceRevCom" runs cutadapt with option -a LPS...RvLPSRevComp 
"""
import os
import subprocess

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'Y':'R', 'R':'Y','S':'S','W':'W','K':'M','M':'K','N':'N','B':'V','V':'B','D':'H','H':'D'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)


def reverse_complement(s):
    return complement(s[::-1])

# List files
fq_files = [f for f in os.listdir(snakemake.params[0]) if f.endswith("_1.fastq.gz")]
if not os.path.exists(snakemake.params[0]+"/reads_discarded_primer/"):
    os.makedirs(snakemake.params[0]+"/reads_discarded_primer/")
if not os.path.exists(snakemake.params[0]+"/primer_removed/"):
    os.makedirs(snakemake.params[0]+"/primer_removed/")
for fw in fq_files:
    sample=fw.replace("_1."+snakemake.params[2],"")
    fw_fq= snakemake.params[0]+"/"+fw
    rv=fw.replace("_1."+snakemake.params[2],"_2."+snakemake.params[2])
    rv_fq= snakemake.params[0]+"/"+rv
    discard_untrimmed=""
    extra_params=snakemake.params[1] 
    if snakemake.params[3] == "PE":
        if "--discard-untrimmed" in snakemake.params[1]:
            discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz --untrimmed-paired-output  "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_2.fastq.gz"
            extra_params=snakemake.params[1].replace("--discard-untrimmed","")
        #print("cutadapt -g "+ primer_by_sample[sample][0] + " -G " + primer_by_sample[sample][1] + " " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log")
        subprocess.run(["cutadapt -g "+ snakemake.config["demultiplexing"]["primers"]["fw_primer"]  + " -G " + snakemake.config["demultiplexing"]["primers"]["rv_primer"]  + " " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz -p "+snakemake.params[0]+"/primer_removed/"+sample+"_2.fastq.gz "+discard_untrimmed +" "+ fw_fq + " " +  rv_fq + " >> "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)
        subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)
        #subprocess.run( ["cutadapt "+ primer_set +" "+snakemake.params[0]+" -o "+snakemake.output[0] + " " + snakemake.input[0]+ ">"+ snakemake.output[1]],stdout=subprocess.PIPE, shell=True)
    elif snakemake.params[3] == "SE":
        if "--discard-untrimmed" in snakemake.params[0]:
            discard_untrimmed=" --untrimmed-output "+snakemake.params[0]+"/reads_discarded_primer/"+sample+"_1.fastq.gz"
            extra_params=snakemake.params[1].replace("--discard-untrimmed","") 
        subprocess.run(["cutadapt -g "+ snakemake.config["demultiplexing"]["primers"]["fw_primer"] +" " +extra_params+" -O "+ snakemake.config["demultiplexing"]["primers"]["min_overlap"] +" -o "+snakemake.params[0]+"/primer_removed/"+sample+"_1.fastq.gz "+ discard_untrimmed + " " + fw_fq + " >> "+ snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log"],stdout=subprocess.PIPE, shell=True)  
        subprocess.run(["grep \"(passing filters)\" "+snakemake.params[0]+"/primer_removed/"+sample+".cutadapt.log | awk '{print \""+sample+"\t\"$5\"\t\"$6}' >> "+snakemake.output[0]],stdout=subprocess.PIPE, shell=True)

