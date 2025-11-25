import os
import subprocess
import shutil
import gzip
import argparse
from argparse import RawTextHelpFormatter
from pathlib import Path
import glob
import subprocess

parser = argparse.ArgumentParser(description="""
This Script validates validates the read assembly process.
It count all the *log files on the (argument) directory
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

parser.add_argument('fastq_directory', metavar='fastq_directory',  type=str, help='Directory containing fastq files from PEAR')
parser.add_argument('-o','--out', dest="out", required = False, default='fastqc_validation.out', type=str, help='Output file with demultiplexing results per sample.')
parser.add_argument('-l','--log', dest="log", required = False, default='sequali_validation.log', type=str, help='Log file of this script.')
parser.add_argument('-q','--qc', dest="qc", choices=['both', 'fastqc', 'sequali'], required = False, default='both', type=str, help='Quality control tool.')
parser.add_argument('-fo','--fastqc-out', dest="fastqc_out", required = False, default='fastqc', type=str, help='Directory for FastQC output files.')
parser.add_argument('-so','--sequali-out', dest="sequali_out", required = False, default='sequali', type=str, help='Directory for Sequali output files.')
parser.add_argument('-cf','--cpus-fastqc', dest="cpus_fastqc", required = False, default=1, type=int, help='Number of threads to run FastQC.')
parser.add_argument('-cs','--cpus-sequali', dest="cpus_sequali", required = False, default=1, type=int, help='Number of threads to run Sequali.')

#FASTQC
#"{config[fastQC][command]} {input} --extract {config[fastQC][extra_params]} -t {config[fastQC][threads]} -o {wildcards.PROJECT}/samples/{wildcards.sample}/qc/"
#SEQUALI
#"sequali --outdir {params.outdir} --html  sequali.html --json sequali.json -t {config[sequali][threads]}  {config[sequali][extra_params]}  {input}"
args = parser.parse_args()
Path(args.fastqc_out).mkdir(parents=True, exist_ok=True)
Path(args.sequali_out).mkdir(parents=True, exist_ok=True)

#sample_log = "#Sample\tReads\tAssembled\tDiscarded\t Not_Assembled\n"
with open(args.out, 'w') as out:
    #out.write("#Sample\tReads\tAssembled\tDiscarded\t Not_Assembled\n")
    #for file in os.listdir(args.fastq_directory):
    for f in Path(args.fastq_directory).glob("*.assembled.fastq"):    
        file=str(f)
        print(file)
        sample=(file.split("/")[-1])[:-16]
        print(sample)
        #sample=str(file)
        #out.write(""+sample)
        if args.qc.lower()=="fastqc":
            subprocess.run( ["fastqc " + file + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["fastqc " + args.fastq_directory+"/"+sample+".unassembled.forward.fastq" + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["fastqc " + args.fastq_directory+"/"+sample+".unassembled.reverse.fastq" + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)  
           
        elif args.qc.lower()=="sequali":
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +file  ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +args.fastq_directory+"/"+sample+".unassembled.forward.fastq"  ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +args.fastq_directory+"/"+sample+".unassembled.reverse.fastq"  ],stdout=subprocess.PIPE, shell=True)  
        else: #both 
            subprocess.run( ["fastqc " + file + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["fastqc " + args.fastq_directory+"/"+sample+".unassembled.forward.fastq" + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["fastqc " + args.fastq_directory+"/"+sample+".unassembled.reverse.fastq" + " --extract -t " + str(args.cpus_fastqc) + " -o " + args.fastqc_out ],stdout=subprocess.PIPE, shell=True)  
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +file  ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +args.fastq_directory+"/"+sample+".unassembled.forward.fastq"  ],stdout=subprocess.PIPE, shell=True)
            subprocess.run( ["sequali --html  sequali."+sample+".html --json sequali."+sample+".json --outdir " + args.sequali_out + " -t "  + str(args.cpus_sequali) +  " " +args.fastq_directory+"/"+sample+".unassembled.reverse.fastq"  ],stdout=subprocess.PIPE, shell=True)  
       
    with open(args.log, "w") as tmplog: 
        tmplog.write("***Read assembly validation***\n") 
 
