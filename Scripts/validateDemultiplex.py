##!/usr/bin/python
# Alejandro Abdala 
# Last recorded update: 07-2025
#
#This Script validates the demultiplexing process within CASCABEL 
#It count all the *_1.fastq.gz files on a directory 


import os
import subprocess
import shutil
import gzip
import argparse
from argparse import RawTextHelpFormatter
from sys import stdin

parser = argparse.ArgumentParser(description="""
This Script validates validates the demultiplexing process within CASCABEL.
It count all the *_1.fastq.gz files on the (argument) directory
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

#parser.add_argument('-m','--mapping-file', metavar='file', required = True, type=str, help='Mapping file to be validated.')
parser.add_argument('fastqc_directory', metavar='fastqc_directory',  type=str, help='Directory containing demultiplexed fastqc files')
parser.add_argument('-r','--raw-read', dest="raw_data", required = True, type=str, help='Raw file used for the demultiplexing process (it could be R1 or R2).')
#parser.add_argument('--is-raw-', required = False,dest="interactive", action=argparse.BooleanOptionalAction,default=True, help='Whether if the script behaviour is interactive or not')
#parser.add_argument('-o','--out', dest="out", required = False, default='summary_demux.txt', type=str, help='Output file with demultiplexing results per sample')
parser.add_argument('--interactive', required = False,dest="interactive", action=argparse.BooleanOptionalAction,default=False, help='Whether if the script behaviour is interactive or not')
parser.add_argument('--tagjump', required = False,dest="tagjump", action=argparse.BooleanOptionalAction,default=False, help='Whether if Cascabel run the tag jump analysis or not.')
parser.add_argument('-o','--out', dest="out", required = False, default='summary_demux.out', type=str, help='Output file with demultiplexing results per sample')
parser.add_argument('-l','--log', dest="log", required = False, default='demux_validation.log', type=str, help='Log file of this script')


args = parser.parse_args()


def count_fastq_gz(file):
    with gzip.open(file, 'rb') as f:
        i=0
        for i, l in enumerate(f):
            pass
    #if i:
    return i + 1 if i else 0



def countReads(fastaFile, isFastq, isGzip):
    command="zcat " if isGzip else "cat "
    if isFastq:
        reads = subprocess.run( [command + fastaFile + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    else:
        reads = subprocess.run( [command + fastaFile + " | grep '^>' | wc -l"],stdout=subprocess.PIPE, shell=True)
    
    totReads =  reads.stdout.decode('utf-8').strip()
    try:
        total = int(totReads)
        if isFastq:
            total = total/4
    except ValueError:
        print("Error counting entries: " + fastaFile)
        total = -1;

    return total

#rd = str(args.raw_data)
#print(rd)
#if args.raw_data.upper().endswith("GZ"): 
#    totRead=countReads(args.raw_data, True, True)
#else:
#    totRead=countReads(args.raw_data, True, False)

t=0
results = []
demux_reads=0
with open(args.out, 'w') as out:
    out.write("Sample\tReads\n")
    for file in os.listdir(args.fastqc_directory):
        if file.endswith("_1.fastq.gz"):
    	    #print(args.fastqc_directory+"/"+file)
    	    sample=file[:-11]
    	    #c=count_fastq_gz(args.fastqc_directory+"/"+file)
    	    c=countReads(args.fastqc_directory+"/"+file, True, True)
    	    demux_reads=c+demux_reads
    	    out.write(sample+"\t"+str(c)+"\n")
    	    results.append(sample+"\t"+str(c))
    	    #print(sample+"\t"+str(c))
    
print("Total: "+str(demux_reads))

if args.raw_data.upper().endswith("GZ"): 
    totReads=countReads(args.raw_data, True, True)
else:
    totReads=countReads(args.raw_data, True, False)

prc = ((totReads - (demux_reads))/totReads)*100

show_menu = True
if  not args.interactive:
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" +"We suggest to review the complete split logs at: "+ args.out+ "\033[0m")
    with open(args.log, "w") as tmplog:
        tmplog.write("Interactive mode off.\n")
        tmplog.write("Total number of input sequences: " + str(totReads) + "\n")
        tmplog.write("Sequences with barcodes in mapping file: " + str(demux_reads) + "\n")
        tmplog.write("Sequences with barcodes not in mapping file: " + str(totReads - (demux_reads)) + " ({0:.2f}".format(prc) +"%) \n")
        tmplog.write("Sequences with those barcodes are dismissed \n")

else:
    cutadapt_dir = "cutadapt".join(args.fastqc_directory.rsplit("demultiplexed", 1))
    while show_menu:
        print("\033[92m" + "Total number of input sequences: " + str(totReads) + "\033[0m")
        print("\033[92m" + "Sequences with barcodes in mapping file: " + str(demux_reads) + "\033[0m")
        print("\033[91m" + "Sequences with barcodes not in mapping file: " + str(totReads - (demux_reads)) + " ({0:.2f}".format(prc) +"%) \033[0m")
        print("\033[91m" + "Sequences with those barcodes are dismiss \033[0m")
        print("\033[92mPlease take a look into complete log files at: \033[93m "+ args.out + " \033[0m")
        print("\033[92mUnassigned reads can be found at file: \033[93m "+ args.fastqc_directory+"/unassigned_R#.fastq.gz \033[0m")
        if args.tagjump:
          mtx = "barcodes".join(args.fastqc_directory.rsplit("demultiplexed", 1))
          print("\033[92mIf too many reads are being discarded, you can have a look at the tag jump analysis result: \033[93m "+ mtx+" \033[0m")
        else:
          print("\033[92mIf too many reads are being discarded, you can set to true the tag jump analysis flag:\"create_tag_jump_analysis\" \033[0m")    
        print("\033[92mYou can press '1' to display detailed sample read assignation\033[0m")


        print("\033[93mDo you want to continue y/n? \033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = " ".join(user_input.split())
        #user_input = user_input[:-1]
        if user_input.upper() == "Y" or user_input.upper() == "YES":
            print("\033[92m" +"The flow goes on!"+ "\033[0m")
            with open(args.log, "w") as tmplog:
                tmplog.write("Workflow continue...")
            exit(0)
        elif user_input.upper() == "N" or user_input.upper() == "NO":
            print("\033[91m" + "Aborting workflow...\033[0m")
            print("\033[92m" + "You can choose to keep or remove current demultiplexed samples. "+"\033[0m")
            print("\033[92m" + "If you remove the results, adjust parameters and restart CASCABEL in order to redo the demultiplexing."+ "\033[0m")
            print("\033[92m" + "If you keep current demultiplexed samples and want to redo it later, you can"+  "\033[0m")
            print("\033[92m" + "restart CASCABEL with the option \"--forcerun cutadapt_demux\" in order to overwrite previous results.  "+ "\033[0m")
            print("\033[93m" + "Do you want to "+ "\033[91m "+"REMOVE" +  "\033[93m "+ "current demultiplexed files?"+ "\033[0m")
            print("\033[92mThis action will remove following directories: \033[93m "+ args.fastqc_directory + " \033[0m")
            print("\033[92mand: \033[93m "+ cutadapt_dir + " \033[0m")
            print("\033[92m" + "Remove directories y/n? "+"\033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = " ".join(user_input.split())
            if user_input.upper() == "Y" or user_input.upper() == "YES":
                print("Cleaning files...")
                shutil.rmtree(args.fastqc_directory)
                shutil.rmtree(cutadapt_dir)

            exit(1)
        elif user_input == "1":
            for sample in results:
                print(sample)
