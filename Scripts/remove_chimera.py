import os
from sys import stdin
import subprocess

try:
    treads = subprocess.run( ["grep '^>' " + snakemake.input[0] + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    totalReads =  treads.stdout.decode('utf-8').strip()
    creads = subprocess.run( ["cat " + snakemake.input[1] + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    chimericReads =  creads.stdout.decode('utf-8').strip()
    prc = (float(chimericReads)/float(totalReads))*100
    print("\033[91m This step can remove possible chimeric sequences \033[0m")
    print("\033[93m Total number of reads: " + totalReads + " \033[0m")
    print("\033[93m Total number of possible chimeras: " + chimericReads + " ({0:.2f}".format(prc) + "%) \033[0m")
    print("\033[92m Do you want to remove chimeric sequences?(y/n): \033[0m")
    if snakemake.config["interactive"] != "F":
        user_input = stdin.readline() #READS A LINE
        user_input = user_input[:-1]
        filter_log = "Total number of possible chimeras: " + chimericReads + " ({0:.2f}".format(prc) + ")%\n\n"
        if user_input.upper() == "Y" or user_input.upper() == "YES":
            subprocess.run( ["filter_fasta.py -f " + snakemake.input[0] + " -s "+ snakemake.input[1] + " -n -o " + snakemake.output[0]], stdout=subprocess.PIPE, shell=True)
            filter_log += "The chimeric sequences were removed with the following command:\n\n"
            filter_log += ":commd:`filter_fasta.py -f " + snakemake.input[0] + " -s "+ snakemake.input[1] + " -n -o " + snakemake.output[0]+"`\n\n"
        else:
            subprocess.run( ["mv " + snakemake.input[0] + " " + snakemake.output[0]], stdout=subprocess.PIPE, shell=True)
            filter_log += "The user didn't remove the chimeric sequences\n\n"
        with open(snakemake.output[1], "w") as out:
            out.write(filter_log)
            out.close()
    else:
        print("\033[93m" +" Interactive mode off \033[0m")
        print("\033[93m" +" Removing chimeras...\033[0m")
        subprocess.run( [snakemake.config["qiime"]["path"]+"filter_fasta.py -f " + snakemake.input[0] + " -s "+ snakemake.input[1] + " -n -o " + snakemake.output[0]], stdout=subprocess.PIPE, shell=True)
        filter_log = "Total number of possible chimeras: " + chimericReads + " ({0:.2f}".format(prc) + ")%\n\n"
        filter_log += "The chimeric sequences were removed with the following command:\n\n"
        filter_log += ":commd:`filter_fasta.py -f " + snakemake.input[0] + " -s "+ snakemake.input[1] + " -n -o " + snakemake.output[0]+"`\n\n"
        with open(snakemake.output[1], "w") as out:
            out.write("Interactive mode off. Automatic chimera removing...\n")
            out.write(str(filter_log))
            out.close()

except Exception as e:
    print("Problem executing script.\nMessage: " + str(e))
