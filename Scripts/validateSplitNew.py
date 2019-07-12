import os
import subprocess
from sys import stdin
import shutil

treads = subprocess.run( ["cat " + snakemake.input.allreads + " | wc -l"],stdout=subprocess.PIPE, shell=True)
totalReads =  treads.stdout.decode('utf-8').strip()
totReads = int(totalReads)/4

sreads = subprocess.run( ["grep '^>' " + snakemake.input.split + " | wc -l"],stdout=subprocess.PIPE, shell=True)
splitReads =  sreads.stdout.decode('utf-8').strip()
spReads = int(splitReads)

sRCreads = subprocess.run( ["grep '^>' " + snakemake.input.splitRC + " | wc -l"],stdout=subprocess.PIPE, shell=True)
splitRCReads =  sRCreads.stdout.decode('utf-8').strip()
sprcReads = int(splitRCReads)
prc = ((totReads - (spReads + sprcReads))/totReads)*100

if snakemake.config["interactive"] == "F":
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" +"We suggest to review the full split logs at: "+ snakemake.input.logSplit+ "\033[0m")
    print("\033[93m" +"As well as: "+ snakemake.input.logSplit+ "\033[0m")
    print("\033[93m" +"And: "+ snakemake.output[0]+ "\033[0m")
    with open(snakemake.output[0], "w") as tmplog:
        tmplog.write("Interactive mode off.\n")
        tmplog.write("Total number of input sequences: " + str(totReads) + "\n")
        tmplog.write("Sequences with barcodes in mapping file: " + str(spReads + sprcReads) + "\n")
        tmplog.write("Sequences with barcodes not in mapping file: " + str(totReads - (spReads + sprcReads)) + " ({0:.2f}".format(prc) +"%) \n")
        tmplog.write("Sequences with those barcodes are dismiss \n")
        tmplog.close()
else:
    print("\033[92m" + "Total number of input sequences: " + str(totReads) + "\033[0m")
    print("\033[92m" + "Sequences with barcodes in mapping file: " + str(spReads + sprcReads) + "\033[0m")
    print("\033[91m" + "Sequences with barcodes not in mapping file: " + str(totReads - (spReads + sprcReads)) + " ({0:.2f}".format(prc) +"%) \033[0m")
    print("\033[91m" + "Sequences with those barcodes are dismiss \033[0m")
    print("\033[92mPlease take a look into complete log files at: "+ snakemake.input.logSplit + " \033[0m")
    print("\033[92mAnd : "+ snakemake.input.logSplitRC + " \033[0m")
    #print("\033[92mYou can run again this step with the flag \"--retain_unassigned_reads\" in order to retain sequences with uncertain barcodes \033[0m")
    print("\033[93mDo you want to continue y/n? \033[0m")
    user_input = stdin.readline() #READS A LINE
    user_input = user_input[:-1]
    if user_input.upper() == "Y" or user_input.upper() == "YES":
        print("\033[92m" +"The flow goes on!"+ "\033[0m")
        with open(snakemake.output[0], "w") as tmplog:
            tmplog.write("Split warning dismissed, user continue...")
            tmplog.close()
    else:
        print("Aborting workflow...")
        print("Cleaning files...")
        shutil.rmtree(snakemake.params[0])
        shutil.rmtree(snakemake.params[1])
        exit(1)
