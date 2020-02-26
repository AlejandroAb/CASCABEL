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
    print("\033[92mPlease take a look into complete log files at: \033[93m "+ snakemake.input.logSplit + " \033[0m")
    print("\033[92mAnd : \033[93m "+ snakemake.input.logSplitRC + " \033[0m")
    print("\033[92mUnassigned reads can be found at file: \033[93m "+ snakemake.input.unassigned + " \033[0m")

    print("\033[93mDo you want to continue y/n? \033[0m")
    user_input = stdin.readline() #READS A LINE
    user_input = " ".join(user_input.split())
    #user_input = user_input[:-1]
    if user_input.upper() == "Y" or user_input.upper() == "YES":
        print("\033[92m" +"The flow goes on!"+ "\033[0m")
        with open(snakemake.output[0], "w") as tmplog:
            tmplog.write("Split warning dismissed, user continue...")
            tmplog.close()
        exit(0)
    else:
        print("\033[91m" + "Aborting workfloW...\033[0m")
        print("\033[92m" + "You can choose to keep or remove current demultiplexed samples. "+"\033[0m")
        print("\033[92m" + "If you remove it, adjust parameters and restart CASCABEL in order to redo the demultiplexing."+ "\033[0m")
        print("\033[92m" + "If you keep current demultiplexed samples and want to redo it later, you can"+  "\033[0m")
        print("\033[92m" + "restart CASCABEL with the option \"--forcerun extract_barcodes\" in order to overwrite previous results.  "+ "\033[0m")
        print("\033[93m" + "Do you want to "+ "\033[91m "+"REMOVE" +  "\033[93m "+ "current demultiplexed files y/n?"+ "\033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = " ".join(user_input.split())
        if user_input.upper() == "Y" or user_input.upper() == "YES":
            print("Cleaning files...")
            shutil.rmtree(snakemake.params[0])
            shutil.rmtree(snakemake.params[1])
        exit(1)
