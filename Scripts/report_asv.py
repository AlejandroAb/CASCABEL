import subprocess
from snakemake.utils import report
import benchmark_utils
from benchmark_utils import readBenchmark
from benchmark_utils import countFasta
from benchmark_utils import countFastaGZ
from benchmark_utils import readSampleDist
from benchmark_utils import make_table
from countData import parseCounts
from seqsChart import createChart

#Parse the total number of counts
#countTxt = parseCounts(snakemake.input.counts)

################################################################################
#                         TOOLS VERSION SECTION                          #
################################################################################
#--fastq
fqv = subprocess.run([snakemake.config["fastQC"]["command"], '--version'], stdout=subprocess.PIPE)
fqVersion = "**" + fqv.stdout.decode('utf-8').strip() + "**"

if snakemake.config["demultiplexing"]["demultiplex"] !=  "F":
   #--qiime extract_barcodes
   ebv = subprocess.run([snakemake.config["qiime"]["path"]+'extract_barcodes.py', '--version'], stdout=subprocess.PIPE)
   ebVersion = ebv.stdout.decode('utf-8')
   ebVersion = "**" + ebVersion[ebVersion.find(":")+1:].strip() + "**"
   #--qiime split_libraries
   spVersion = "**TBD**"
   spv = subprocess.run([snakemake.config["qiime"]["path"]+'split_libraries_fastq.py', '--version'], stdout=subprocess.PIPE)
   spVersion = spv.stdout.decode('utf-8')
   if "Version" in spVersion:
       spVersion = "**" + spVersion[spVersion.find(":")+1:].strip() + "**"
else:
   ebVersion = "**NA**"
   SPvERSION = "**NA**"

vsearchVersion = "**TBD**"
vsearchV = subprocess.run(['vsearch', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
vsearchVersion = "**" + vsearchV.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"


#--qiime identify_chimeric_seqs
icVersion = "**TBD**"
icv = subprocess.run([snakemake.config["qiime"]["path"]+'identify_chimeric_seqs.py', '--version'], stdout=subprocess.PIPE)
icVersion = icv.stdout.decode('utf-8')
if "Version" in icVersion:
    icVersion = "**" + icVersion[icVersion.find(":")+1:].strip() + "**"
#--pear
try:
    pearv = subprocess.run( [snakemake.config["pear"]["command"]+" -h | grep 'PEAR v'"], stdout=subprocess.PIPE, shell=True)
    pearversion = "**" + pearv.stdout.decode('utf-8').strip() + "**"
except Exception as e:
    pearversion = "Problem reading version"

#--cutadapt
if snakemake.config["cutAdapters"] == "T":
    #cutv = subprocess.run(['/export/data/aabdala/.local/bin/cutadapt', '--version'], stdout=subprocess.PIPE)
    #cutVersion = "**cutadapt v" + cutv.stdout.decode('utf-8').strip() + "**"
    cutVersion = "cutadapt v TBD"

################################################################################
#                          Chimera check                                       #
################################################################################
removeChimeras = False
if snakemake.config["chimera"]["search"] == "T":
    ################################################################################
    #                       Read log file from remove_chimera.py                   #
    # After search for chimera, user have the option to remove them or not. If the #
    # user decides to remove the chimera, the executed command is stored on the log#
    # file, otherwise it stores a message indicating the user decision.            #
    ################################################################################
    chimera_log = ""
    try:
        with open(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/chimera/chimera.log") as chimlog:
            for line in chimlog:
                chimera_log += line
            chimlog.close()
    except FileNotFoundError:
        chiemra_log = "No Log for identify_chimeric_seqs.py"
    if "The chimeric sequences were removed" in chimera_log:
        removeChimeras = True


################################################################################
#                         Benchmark Section                                    #
# This section is to generate a pre-formatted text with the benchmark info for #
# All the executed rules.                                                      #
################################################################################
fqBench = readBenchmark(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/qc/fq.benchmark")
pearBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/pear.benchmark")
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    barBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.benchmark")
    splitLibsBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/splitLibs.benchmark")
    splitLibsRCBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/splitLibs.benchmark")
    combineBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/combine_seqs_fw_rev.benchmark")
else:
    combineBench=pearBench #THIS IS ONLY FOR TESTING REMOVE!!! 
rmShorLongBench = ""
if snakemake.config["ANALYSIS_TYPE"] == "OTU": 
    rmShorLongBench = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/filter.benchmark")
demultiplexFQBench=""
if snakemake.config["demultiplexing"]["demultiplex"] == "T" and snakemake.config["demultiplexing"]["create_fastq_files"] == "T":
    demultiplexFQBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/demultiplex_fq.benchmark")

################################################################################
#                           Compute Counts                                     #
################################################################################
if snakemake.config["gzip_input"] == "F":
    rawCounts = countFasta(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/rawdata/fw.fastq", True);
else:
    rawCounts = countFastaGZ(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/rawdata/fw.fastq.gz", True);
rawCountsStr= str(int(rawCounts))
#rawCountsStr= '{0:g}'.format(float(rawCounts))
#-peared
pearedCounts = 0
if snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T":
    pearedCounts = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq", True);
else:
    pearedCounts = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.UNPAIRED.fastq", True);

pearedCountsStr=str(int(pearedCounts))
#pearedCountsStr='{0:g}'.format(float(pearedCounts))
prcPeared = "{:.2f}".format(float((pearedCounts/rawCounts)*100))
#-dumultiplex
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    fwAssignedCounts = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.no_unassigneds.fna", False)
    rvAssignedCounts = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/seqs.no_unassigneds.fna", False)
    prcFwAssigned = "{:.2f}".format(float((fwAssignedCounts/pearedCounts)*100))
    prcRvAssigned = "{:.2f}".format(float((rvAssignedCounts/pearedCounts)*100))
    totalAssigned = fwAssignedCounts + rvAssignedCounts
    prcPearedAssigned = "{:.2f}".format(float((totalAssigned/pearedCounts)*100))
    prcRawAssigned = "{:.2f}".format(float((totalAssigned/rawCounts)*100))
else: 
    totalAssigned = pearedCounts
    prcPearedAssigned = "{:.2f}".format(float((totalAssigned/pearedCounts)*100))
    prcRawAssigned = "{:.2f}".format(float((totalAssigned/rawCounts)*100))

#--cutadapt
cutSequences = False
if snakemake.config["cutAdapters"] == "T":
    sequenceNoAdapters = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.fna", False)
    if (totalAssigned - sequenceNoAdapters) > 0:
        cutSequences = True
        prcCut = "{:.2f}".format(float((sequenceNoAdapters/totalAssigned)*100))
        prcCutRaw = "{:.2f}".format(float((sequenceNoAdapters/rawCounts)*100))

if removeChimeras:
    sequenceNoChimeras = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered_nc.fasta", False)
    prcChim = "{:.2f}".format(float((sequenceNoChimeras/totalAssigned)*100))
    prcChimRaw = "{:.2f}".format(float((sequenceNoChimeras/rawCounts)*100))
    if cutSequences:
        prcChimCut = "{:.2f}".format(float((sequenceNoChimeras/sequenceNoAdapters)*100))
#out="{PROJECT}/runs/{run}/{sample}_data/"
trimmedCounts=1;
if snakemake.config["ANALYSIS_TYPE"] == "OTU":
    trimmedCounts = countFasta(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.fasta", False)
    prcTrimmedSplit ="{:.2f}".format(float((trimmedCounts/totalAssigned)*100))
    prcTrimmedRaw= "{:.2f}".format(float((trimmedCounts/rawCounts)*100))
if cutSequences:
    prcTrimmedCut="{:.2f}".format(float((trimmedCounts/sequenceNoAdapters)*100))
#if removeChimeras:
#    prcTrimmedChimera="{:.2f}".format(float((trimmedCounts/sequenceNoChimeras)*100))
try:
    samplesLib = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.dist.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    samplesLibInt = int(samplesLib.stdout.decode('utf-8').strip())
except Exception as e:
    totalReads = "Problem reading outputfile"
################################################################################
#                         Generate sequence amounts chart                      #
################################################################################
numbers=[rawCounts,pearedCounts];
labels=["Raw", "Assembled"];
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    numbers.append(totalAssigned)
    labels.append("Demultiplexed")
if snakemake.config["cutAdapters"] == "T":
    numbers.append(sequenceNoAdapters)
    labels.append("Cutadapt")
#if snakemake.config["ANALYSIS_TYPE"] == "OTU":
#    numbers.append(trimmedCounts)
#    labels.append("Length filtering")
#if snakemake.config["chimera"]["search"] == "T":
#    numbers.append(sequenceNoChimeras)
#    labels.append("No Chimera")
createChart(numbers, tuple(labels),snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files/sequence_numbers."+snakemake.wildcards.sample+".png")
################################################################################
#                          Chimera check                                       #
################################################################################
variable_refs=""
#if snakemake.config["chimera"]["search"] == "T" and snakemake.config["chimera"]["method"] == "usearch61":
#    variable_refs+= ".. [usearch61] Edgar RC. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-2461.\n\n"
#else: 
#    variable_refs+= ".. [uchime] Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics, 27 (16): 2194-2200. doi:10.1093/bioinformatics/btr381. \n\n"
quimeraStr = ""
#if snakemake.config["chimera"]["search"] == "T":
#    quimeraStr="Identify Chimera\n-------------------\n\n"
#    quimeraStr+="Identify possible chimeric sequences (sequences generated due to the PCR amplification of multiple templates or parent sequences).\n\n"
#    if snakemake.config["chimera"]["method"] == "usearch61":
#        quimeraStr += ":red:`Tool:` [QIIME]_ - identify_chimeric_seqs.py\n\n"
#        quimeraStr += ":red:`Version:` "+ icVersion +"\n\n"
#        quimeraStr += ":red:`Method:` [usearch61]_ \n\n"
#    else:
#        quimeraStr += ":red:`Tool:` [Vsearch]_ - vsearch\n\n"
#        quimeraStr += ":red:`Version:` "+ vsearchVersion +"\n\n"        
#        quimeraStr += ":red:`Method:` "+ str(snakemake.config["chimera"]["method"]) +" - uses [uchime]_ \n\n"    
#    quimeraStr += "**Command:**\n\n"
#    if snakemake.config["chimera"]["method"] == "usearch61":
#         quimeraStr+=":commd:`identify_chimeric_seqs.py -m "+ str(snakemake.config["chimera"]["method"])+" -i "+ str(snakemake.wildcards.PROJECT)+"/runs/"+str(snakemake.wildcards.run)+"/"+str(snakemake.wildcards.sample)+"_data/seqs_fw_rev_accepted.fna "+str(snakemake.config["chimera"]["extra_params"])
#         quimeraStr+=" -o "+ str(snakemake.wildcards.PROJECT)+"/runs/"+str(snakemake.wildcards.run)+"/"+str(snakemake.wildcards.sample)+"_data/chimera/` \n\n"
#    else:
#         quimeraStr+=":commd:`vsearch --"+ str(snakemake.config["chimera"]["method"])+" "+ str(snakemake.wildcards.PROJECT)+"/runs/"+str(snakemake.wildcards.run)+"/"+str(snakemake.wildcards.sample)+"_data/seqs_fw_rev_accepted.fna --threads "+ str(snakemake.config["chimera"]["threads"]) +" " +str(snakemake.config["chimera"]["extra_params"])   
#         quimeraStr+=" --uchimeout "+ str(snakemake.wildcards.PROJECT)+"/runs/"+str(snakemake.wildcards.run)+"/"+str(snakemake.wildcards.sample)+"_data/chimera/chimeras.summary.txt` \n\n"
#    quimeraStr+="**Output files:**\n\n"
#    if snakemake.config["chimera"]["method"] == "usearch61":
#        quimeraStr+=":green:`- File with the possible chimeric sequences:` "+str(snakemake.wildcards.PROJECT)+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/chimera/chimeras.txt\n\n"
#    else:
#        quimeraStr+=":green:`- File with the possible chimeric sequences:` "+str(snakemake.wildcards.PROJECT)+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/chimera/chimeras.summary.txt\n\n"
#    identifyChimeraBench=readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/chimera/chimera.benchmark")
#    quimeraStr+=identifyChimeraBench
#    quimeraStr+=chimera_log
#    if removeChimeras:
#        quimeraStr+=":red:`Reads after remove chimeric sequences:` "+ str(sequenceNoChimeras)+"\n\n"
#        quimeraStr+=":red:`Percentage of reads vs raw reads:` "+ str(prcChimRaw) + "%\n\n"
#        quimeraStr+=":red:`Percentage of reads vs demultiplexed reads:` "+ str(prcChim) + "%\n\n"
#        if cutSequences:
#            quimeraStr+=":red:`Percentage of reads vs cutadapt:` "+ str(prcChimRaw) + "%\n\n"




################################################################################
#                           Peared FastQC                                     #
################################################################################
fastQCPearStr = ""
if snakemake.config["fastQCPear"] == "T":
    fastQCPearStr = "Peared FastQC Analysis\n------------------------\n\n" # title
    fastQCPearStr += "Check the quality of the reads after assembly.\n\n"
    fastQCPearStr += ":red:`Tool:` [FastQC]_\n\n"
    fastQCPearStr += ":red:`Version:` "+ fqVersion +"\n\n"
    fastQCPearStr += "**Command:**\n\n"
    fastQCPearStr += ":commd:`fastqc "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq --extract -o  "+ snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/qc`\n\n"
    fastQCPearStr += "**Output files:**\n\n:green:`- FastQC report:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/qc/seqs.assembled_fastqc.html FQ_Report_ \n\n"
    fastQCPearStr += ".. _FQ_Report: peared/qc/seqs.assembled_fastqc.html \n\n"
    fastQCPearStrBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/qc/fq.benchmark")
    fastQCPearStr += fastQCPearStrBench

################################################################################
#                           Extract Barcode                                    #
################################################################################
extractBCStr = ""
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    extractBCStr ="Extract barcodes\n-----------------\n\n"
    extractBCStr +="Extract the barcodes used to identify individual samples.\n\n"
    extractBCStr +=":red:`Tool:` [QIIME]_ - extract_barcodes.py\n\n"
    extractBCStr +=":red:`Version:` "+ebVersion+"\n\n"
    extractBCStr +="**Command:**\n\n"
    extractBCStr +=":commd:`extract_barcodes.py -f "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq -c "+str(snakemake.config["ext_bc"]["c"])+ " " + str(snakemake.config["ext_bc"]["bc_length"])+ " " + snakemake.config["ext_bc"]["extra_params"] + " -o "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/`\n\n"
    extractBCStr +="**Output files:**\n\n"
    extractBCStr +=":green:`- Fastq file with barcodes:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq\n\n"
    extractBCStr +=":green:`- Fastq file with the reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/reads.fastq\n\n"
    extractBCStr +=barBench
################################################################################
#                           CORRECT Barcodeds                                  #
################################################################################
correctBCStr = ""
bcFile="barcodes.fastq"
if snakemake.config["demultiplexing"]["demultiplex"] != "F" and snakemake.config["bc_mismatch"]:
    correctBCStr = "Correct Barcodes\n--------------------\n"
    correctBCStr += "Try to correct the barcode from unassigned reads.\n\n"
    correctBCStr += "Maximum number of mismatches **"  + str(snakemake.config["bc_mismatch"]) + "**.\n\n"
    correctBCStr +=":red:`Tool:` CASCABEL's R script\n\n"
    correctBCStr +="**Command:**\n\n"
    correctBCStr += ":commd:`Rscript Scripts/errorCorrectBarcodes.R $PWD "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq "  + str(snakemake.config["bc_mismatch"]) + "`\n\n"
    correctBCStr += "**Output file:**\n\n:green:`- Barcode corrected file:` "+snakemake.wildcards.PROJECT+ "/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq_corrected\n\n"
    correctBarBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes_corrected.benchmark")
    correctBCStr += correctBarBench
    bcFile="barcodes.fastq_corrected"

splitStr = ""
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    splitStr+="Demultiplexing\n"
    splitStr+="----------------\n"
    splitStr+="For library splitting, also known as demultiplexing, Cascabel performs several steps to assign fragments in the original as well as reverse orientation to the correct sample.\n\n"
    splitStr+="Split samples from Fastq file\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+=":red:`Tool:` [QIIME]_ - split_libraries_fastq.py\n\n"
    splitStr+=":red:`version:` "+ spVersion+"\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`split_libraries_fastq.py -m "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -i "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/reads.fastq -o  "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs -b "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/"+bcFile+" -q "+str(snakemake.config["split"]["q"])+" -r "+str(snakemake.config["split"]["r"])+" --retain_unassigned_reads "+str(snakemake.config["split"]["extra_params"])+" --barcode_type "+str(snakemake.config["split"]["barcode_type"])+"`\n\n"
    splitStr+=splitLibsBench

    splitStr+="Retain assigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`cat "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.fna | grep -P -A1 \"(?!>Unass)^>\" | sed '/^--$/d' > "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.assigned.fna`\n\n"

    splitStr+="Create file with only unassigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`cat "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.fna | grep \"^>Unassigned\" |  sed 's/>Unassigned_[0-9]* /@/g' | sed 's/ .*//' | grep -F -w -A3  -f - "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq |  sed '/^--$/d' >"+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/unassigned.fastq`\n\n"

    splitStr+="Reverse complement unassigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+=":red:`Tool:` [Vsearch]_\n\n"
    splitStr+=":red:`version:`  "+vsearchVersion+"\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`vsearch --fastx_revcomp "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/unassigned.fastq  --fastqout "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/unassigned.reversed.fastq`\n\n"


    splitStr+="Barcode extraction for reverse complemented, unassigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr +=":red:`Tool:` [QIIME]_ - extract_barcodes.py\n\n"
    splitStr +=":red:`Version:` "+ebVersion+"\n\n"
    splitStr+="**Command:**\n\n"
    splitStr +=":commd:`extract_barcodes.py -f "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/unassigned.reversed.fastq -c "+str(snakemake.config["ext_bc"]["c"])+" "+str(snakemake.config["ext_bc"]["bc_length"])+" "+snakemake.config["ext_bc"]["extra_params"]+" -o "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/`\n\n"

    if snakemake.config["bc_mismatch"]:
        splitStr += "Correct reverse complemented barcodes \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        splitStr += "Maximum number of mismatches **"  + str(snakemake.config["bc_mismatch"]) + "**.\n\n"
        splitStr +=":red:`Tool:` CASCABEL's R script\n\n"
        splitStr +="**Command:**\n\n"
        splitStr += ":commd:`Rscript Scripts/errorCorrectBarcodes.R $PWD "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/barcodes.fastq_corrected "  + str(snakemake.config["bc_mismatch"]) + "`\n\n"
        splitStr += "**Output file:**\n\n:green:`- Barcode corrected file:` "+snakemake.wildcards.PROJECT+ "/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes/barcodes.fastq_corrected\n\n"
        splitStrBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/barcodes_corrected.benchmark")
        splitStr += splitStrBench+"\n\n"

    splitStr +="Split reverse complemented reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr +=":red:`Tool:` [QIIME]_ - extract_barcodes.py\n\n"
    splitStr +=":red:`Version:` "+ebVersion+"\n\n"
    splitStr+="**Command:**\n\n"
    splitStr +=":commd:`split_libraries_fastq.py -m "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -i "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/barcodes_unassigned/reads.fastq -o "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC -b "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+str(snakemake.wildcards.sample)+"_data/barcodes_unassigned/"+bcFile+" -q "+str(snakemake.config["split"]["q"])+" -r "+str(snakemake.config["split"]["r"])+" "+str(snakemake.config["split"]["extra_params"])+" --barcode_type "+str(snakemake.config["split"]["barcode_type"])+"`\n\n"
    splitStr +=splitLibsBench+"\n\n"

    splitStr +="**Output files:**\n\n"
   # splitStr +=":green:`- FW reads fasta file with new header:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.assigned.fna\n\n"
    splitStr +=":green:`- Text histogram with the length of the fw reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/histograms.txt\n\n"
    splitStr +=":green:`- Log file for the fw reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/split_library_log.txt\n\n"
   # splitStr +=":green:`- RV reads fasta file with new header:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/seqs.assigned.fna\n\n"
    splitStr +=":green:`- Text histogram with the length of the rv reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/histograms.txt\n\n"
    splitStr +=":green:`- Log file for the rv reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/split_library_log.txt\n\n"
    splitStr +=":green:`- Fasta file with unassigned reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/seqs.unassigned.fna\n\n"
    splitStr +=":red:`Number of reads assigned on FW:` "+str(fwAssignedCounts)+" = "+str(prcFwAssigned)+"% of the peared reads\n\n"
    splitStr +=":red:`Number of reads assigned on RVC:` "+str(rvAssignedCounts)+" = "+str(prcRvAssigned)+"% of the peared reads\n\n"

################################################################################
#                           Single FastQ creation                              #
################################################################################
demultiplexFQ = ""
if snakemake.config["demultiplexing"]["demultiplex"] == "T" and snakemake.config["demultiplexing"]["create_fastq_files"] == "T":
    demultiplexFQ = "Generate single sample fastq files\n------------------------------------------\n\n" # title
    demultiplexFQ += "Create single fastq files per samples (based on the raw data without applying any filtering).\n\n"
    demultiplexFQ +=":red:`Tool:` CASCABEL's Java program\n\n"
    demultiplexFQ += "**Command:**\n\n"
    demultiplexFQ += ":commd:`"+snakemake.config["java"]["command"] + " -cp Scripts DemultiplexQiime --fasta -d "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna -o "+ snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/ "
    ext=".gz"
    if snakemake.config["gzip_input"] == "F":
        ext=""
    demultiplexFQ += "-r1 "+snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/rawdata/fw.fastq"+ext+" -r2 "+snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/rawdata/fw.fastq"+ext+"`\n\n"
    demultiplexFQ += "**The demultiplexed files are located at:**\n\n:green:`- demultiplexed directory:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/\n\n"
    demultiplexFQ += ":green:`- summary file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/demultiplexed/summary.txt\n\n"
    demultiplexFQ += demultiplexFQBench


################################################################################
#                           Combine FW and Reverse reads                       #
################################################################################

combineFR = ""
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    combineFR = "Combine reads\n---------------------------------\n\n" # title
    combineFR += "Concatenate forward and reverse reads.\n\n"
    combineFR += "**Command:**\n\n"
    combineFR += ":commd:`cat "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibs/seqs.assigned.fna "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/splitLibsRC/seqs.assigned.fna > "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna`\n\n"
    combineFR +="**Output files:**\n\n"
    combineFR +=":green:`- Fasta file with combined reads:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna\n\n"
    combineFR +=":red:`- Total number of acepted reads:` " +str(totalAssigned)+ " = "+ str(prcPearedAssigned)+ "% of the peared reads or "+str(prcRawAssigned)+"% of the raw reads.\n\n"
    combineFR += combineBench
    
################################################################################
#                          Cut adapters                                        #
################################################################################
cutAdaptStr = ""
if snakemake.config["cutAdapters"] == "T":
    cutAdaptStr = "Remove sequence primers\n------------------------\n\n" # title
    cutAdaptStr +="Remove the adapters / primers from the reads.\n\n"
    demultiplexFQ +=":red:`Tool:` [Cutadapt]_\n\n"
    cutAdaptStr += ":red:`Version:` "+cutVersion+"\n\n"
    cutAdaptStr += "**Command:**\n\n"
    cutAdaptStr += ":commd:`cutadapt "+ str(snakemake.config["cutadapt"]["adapters"])+" " + str(snakemake.config["cutadapt"]["extra_params"]) + "-o " + snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.fna`\n\n"
    cutAdaptStr +=  snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna > " +  snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.log`\n\n"
    cutAdaptStr += "**Output files:**\n\n:green:`- Reads without adapters:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.fna\n\n"
    if cutSequences:
        cutAdaptStr += ":red:`Total number of reads after cutadapt:` "+ str(sequenceNoAdapters) + " = " + str(prcCut) + "% of the assigned reads or "+ str(prcCutRaw)+"% of the total reads"
    cutAdaptBench =readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/cutadapt.benchmark")
    cutAdaptStr += cutAdaptBench+"\n\n"
################################################################################
#                          Counts for too long too shorts                      #
################################################################################
#trimmedStr =  ":red:`Total number of reads after trimming:` "+str(trimmedCounts)+ "="+ str(prcTrimmedSplit)+"% of the demultiplexed reads or " + str(prcTrimmedRaw) + "% of the raw reads\n\n"
trimmedStr = ""
if snakemake.config["ANALYSIS_TYPE"] == "OTU":
    trimmedStr =  ":red:`Total number of reads after length filtering:` "+str(trimmedCounts)+ "\n\n"
    trimmedStr += ":red:`Percentage of reads vs raw reads:` "+str(prcTrimmedRaw)+"%\n\n"
    trimmedStr+=":red:`Percentage of reads vs demultiplexed reads:` " + str(prcTrimmedSplit) + "%\n\n"
    if cutSequences:
        trimmedStr+=":red:`Percentage of reads after cutadapt:` "+ str(prcTrimmedCut) + "%\n"
#if removeChimeras:
#    trimmedStr+=":red:`Percentage of reads after remove chimeras vs trimmed reads:` "+ str(prcTrimmedChimera) + "%\n"


#bcValidationBench =readBenchmark(snakemake.wildcards.PROJECT+"/metadata/bc_validation/"+snakemake.wildcards.sample+"/validation.benchmark")
################################################################################
#                     Remove too short and too long reads                      #
#  This rule creates a temporary file with the short and long values choosed   #
#  by the user in order to remove the reads. The file filter.log contains the  #
#  minimun expected length for the reads followed by the maximun length tab    #
#  separated (shorts <TAB> longs)                                              #
################################################################################
if snakemake.config["ANALYSIS_TYPE"] == "OTU": 
    shorts = str(snakemake.config["rm_reads"]["shorts"])
    longs = str(snakemake.config["rm_reads"]["longs"])
    with open(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/filter.log") as trimlog:
        for line in trimlog:
            tokens = line.split("\t")
            if len(tokens)>2:
                shorts = tokens[1]
                longs = tokens[2]
################################################################################
#                       FInal Counts                              #
################################################################################
countTxt="Following you can see the final read counts: \n\n"
fileData = []
headers = []
data =[]
headers.append("File description")
headers.append("Location")
headers.append("Number of reads")
headers.append("Prc(%) vs raw")
fileData.append(headers)
#raw
data.append("Raw reads")
data.append(snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample+"/rawdata/\*.fq")
data.append(str(rawCounts))
data.append("{:.2f}".format(float((rawCounts/rawCounts)*100))+"%")
fileData.append(data)
data=[]
#pear
data.append("Assembled reads")
data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/seqs.assembled.fastq")
data.append(str(pearedCounts))
data.append("{:.2f}".format(float((pearedCounts/rawCounts)*100))+"%")
fileData.append(data)
data=[]
#splitted
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    data.append("Demultiplexed reads")
    data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted.fna")
    data.append(str(totalAssigned))
    data.append("{:.2f}".format(float((totalAssigned/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
#adapters
if snakemake.config["cutAdapters"] == "T":
    data.append("Adapter removed")
    data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_accepted_no_adapters.fna")
    data.append(str(sequenceNoAdapters))
    data.append("{:.2f}".format(float((sequenceNoAdapters/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
#length filtered
if snakemake.config["ANALYSIS_TYPE"] == "OTU":
    data.append("Length filtered")
    data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.fasta")
    data.append(str(trimmedCounts))
    data.append("{:.2f}".format(float((trimmedCounts/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
#chimera
if snakemake.config["chimera"]["search"] == "T":
    data.append("Non chimeric reads")
    data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered_nc.fasta")
    data.append(str(sequenceNoChimeras))
    data.append("{:.2f}".format(float((sequenceNoChimeras/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
countTxt += make_table(fileData)
################################################################################
#                       Sample distribution chart                              #
################################################################################

sampleDistChart = ""
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    dist_table = readSampleDist(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.dist.txt",trimmedCounts,samplesLibInt)
    sampleDistChart = "Sample distribution\n--------------------------------------\n\n" # title
    sampleDistChart += dist_table + "\n\n"
    sampleDistChart += ".. image:: report_files/seqs_fw_rev_filtered."+snakemake.wildcards.sample+".dist.png\n\n"
    sampleDistChart +="The previous chart shows the number of clean reads per sample. The bars are sorted from left to right, according to the metadata input file.\n\n"
    sampleDistChart +="**To see more details about the number of reads per sample in this library, please refer to the file:** "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.dist.txt\n\n"


################################################################################
#                       User description section                               #
################################################################################
desc = snakemake.config["description"]
txtDescription = ""
if len(desc) > 0:
    txtDescription = "\n**User description:** "+desc+"\n"

################################################################################
#                                Report                                        #
################################################################################

report("""
Amplicon Analysis Report for Library: {snakemake.wildcards.sample}
=====================================================================
    .. role:: commd
    .. role:: red
    .. role:: green

**CASCABEL** is designed to run amplicon sequence analysis across single or multiple read libraries.

The objective of this pipeline is to create different output files which allow the user to explore data in a simple and meaningful way, as well as facilitate downstream analysis, based on the generated output files.

Another aim of **CASCABEL** is also to encourage the documentation process, by creating this report in order to assure data analysis reproducibility.

{txtDescription}

Following you can see all the steps that were taken in order to get the final results of the pipeline.

Raw Data
---------
The raw data for this library can be found at:

:green:`- FW raw reads:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq

:green:`- RV raw reads:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq

:red:`Number of total reads:` {rawCountsStr}

Quality Control
------------------
Evaluate quality on raw reads.

:red:`Tool:` [FastQC]_

:red:`Version:` {fqVersion}

**Command:**

:commd:`fastqc {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq --extract -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/qc/`

You can follow the links below, in order to see the complete FastQC report:

:green:`- FastQC for sample {snakemake.wildcards.sample}_1:` FQ1_

    .. _FQ1: ../../../samples/{snakemake.wildcards.sample}/qc/fw_fastqc.html

:green:`- FastQC for sample {snakemake.wildcards.sample}_2:` FQ2_

    .. _FQ2: ../../../samples/{snakemake.wildcards.sample}/qc/rv_fastqc.html

{fqBench}


Read pairing
----------------
Align paired end reads and merge them into one single sequence in case they overlap.

:red:`Tool:` [PEAR]_

:red:`version:` {pearversion}

**Command:**

:commd:`pear -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/fw.fastq -r {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.sample}/rawdata/rv.fastq -t {snakemake.config[pear][t]} -v {snakemake.config[pear][v]} -j {snakemake.config[pear][j]} -p {snakemake.config[pear][p]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs > {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs.assembled.fastq`

**Output files:**

:green:`- Merged reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/seqs.assembled.fastq

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/{snakemake.wildcards.sample}_data/peared/pear.log

:red:`Number of peared reads:` {pearedCountsStr} =  {prcPeared}%

{pearBench}

{fastQCPearStr}

{extractBCStr}

{correctBCStr}

{splitStr}

{demultiplexFQ}

{combineFR}

{cutAdaptStr}


{quimeraStr}


{sampleDistChart}


Final counts
-------------

{countTxt}

.. image:: report_files/sequence_numbers.{snakemake.wildcards.sample}.png

ASV report
---------------------------

Cascabel report on downstream analyses in combination with multiple libraries (if supplied), can be found at the following link: asv_report_ ({snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv_report_dada2.html)

    .. _asv_report: asv_report_dada2.html

References
------------------

.. [FastQC] FastQC v0.11.3. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data

.. [PEAR] PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593

.. [QIIME] QIIME. Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Gonzalez Pena A, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. 2010. QIIME allows analysis of high-throughput community sequencing data. Nature Methods 7(5): 335-336.

.. [Cutadapt] Cutadapt v1.15 .Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

.. [Vsearch] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584


{variable_refs}


""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
