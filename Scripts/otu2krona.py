import os
import subprocess
import re

samples = snakemake.config["krona"]["samples"]

tuplesToPrint = []
sampleList = []
if samples.strip() != "all":
    sampleList = [x.strip() for x in re.split(',|;',samples)]

with open(snakemake.input[0]) as otuTxt:
    for line in otuTxt:
        if "#OTU" in line:
            allSamps = line.rstrip('\n').split('\t')
            for index, samp in enumerate(allSamps):
                if index > 0 and index < len(allSamps)-1:
                    if samples.strip() == "all" or samp in sampleList:
                            tuplesToPrint.append((index+1,samp))
            break
    otuTxt.close()
cmmd = snakemake.config["krona"]["command"] + " "
for samp in tuplesToPrint:
    #print("cat "+snakemake.input[0] + " | grep -v \"^#\" | cut -f"+str(samp[0])+","+str(len(allSamps))+" | sed 's/;/\\t/g' | sed 's/*/no_rank/g' > "+samp[1]+".txt")
    subprocess.run(["cat "+snakemake.input[0] + " | grep -v \"^#\" | cut -f"+str(samp[0])+","+str(len(allSamps))+" | grep -v \"^0\" | sed \'s/;/\\t/g\' | sed \'s/*/no_rank/g\' > "+snakemake.params[0]+samp[1]+".krona.txt"],stdout=subprocess.PIPE, shell=True)
    cmmd+=snakemake.params[0]+samp[1]+".krona.txt,"+samp[1] + " "
cmmd+=" -o " + snakemake.output[0] + " -n root " + snakemake.config["krona"]["extra_params"]
out = subprocess.run([cmmd],stdout=subprocess.PIPE, shell=True)
print("Krona report done!")

print("Cleaning intermediate files...")

for samp in tuplesToPrint:
    subprocess.run(["rm -f "+snakemake.params[0]+samp[1]+".krona.txt"],stdout=subprocess.PIPE, shell=True)

exit(0)
