import subprocess
from snakemake.utils import report
import benchmark_utils 
from benchmark_utils import countTxt
from benchmark_utils import readBenchmark
from benchmark_utils import countFasta
from benchmark_utils import countFastaGZ
from benchmark_utils import readSampleDist
from benchmark_utils import make_table
from countData import parseCounts
from seqsChart import createChart

#Parse the total number of counts
#countTxt = parseCounts(snakemake.input.counts)

########################################################
#                CASCABEL version                      #
########################################################
version = ""
try:
  with open('cascabel.version') as f:
    version = f.readline().strip('\n')
except FileNotFoundError:
        print("file version missing: ../cascabel.version\nYou can see README file for Cascabel version.")
        version = "see README file"   

########################################################
#      Base directories for paths                      #
########################################################
base_sample= snakemake.wildcards.PROJECT+"/samples/"+snakemake.wildcards.sample
base_run= snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run
base_sample_data=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data"
base_barcodes=base_sample+"/barcodes"  if snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T" else base_sample+"/barcodes_unpaired"
base_split=base_sample+"/splitLibs"  if snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T" else base_sample+"/splitLibs_unpaired"
base_demultiplexed=base_sample+"/demultiplexed"  if snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T" else base_sample+"/demultiplexed/unpaired"
#snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/peared/pear.benchmark")
################################################################################
#                         TOOLS VERSION SECTION                          #
################################################################################
#--fastq
fqVersion=""
if snakemake.config["QC"].lower() ==  "qc" or snakemake.config["QC"].lower() ==  "both":
    fqv = subprocess.run([snakemake.config["fastQC"]["command"], '--version'], stdout=subprocess.PIPE)
    fqVersion = "**" + fqv.stdout.decode('utf-8').strip() + "**"
sequaliVersion = ""
if snakemake.config["QC"].lower() ==  "sequali" or snakemake.config["QC"].lower() ==  "both":
    sqv = subprocess.run(['sequali', '--version'], stdout=subprocess.PIPE)
    sequaliVersion = "**" + sqv.stdout.decode('utf-8').strip() + "**"


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
cutVersion = "**TBD**"
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg" or snakemake.config["primers"]["remove"].lower() != "f":
    cutv = subprocess.run(['cutadapt', '--version'], stdout=subprocess.PIPE)
    cutVersion = "**cutadapt v" + cutv.stdout.decode('utf-8').strip() + "**"
    #cutVersion = "cutadapt v TBD"

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
        with open(base_sample_data+"/chimera/chimera.log") as chimlog:
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
#fqBench = readBenchmark(base_sample+"/qc/fq.benchmark")
pearBench =readBenchmark(base_sample+"/peared/pear.benchmark")
if snakemake.config["demultiplexing"]["demultiplex"] != "F":
    barBench =readBenchmark(base_barcodes+"/barcodes.benchmark")
    splitLibsBench =readBenchmark(base_split+"/splitLibs.benchmark")
    #splitLibsRCBench =readBenchmark(base_sample_data+"/splitLibsRC/splitLibs.benchmark")
   # combineBench =readBenchmark(base_sample_data+"/combine_seqs_fw_rev.benchmark")
else:
    combineBench=pearBench #THIS IS ONLY FOR TESTING REMOVE!!! 
rmShorLongBench =readBenchmark(base_sample_data+"/filter.benchmark")
demultiplexFQBench=""
if snakemake.config["demultiplexing"]["demultiplex"] == "T" and snakemake.config["demultiplexing"]["create_fastq_files"] == "T":
    demultiplexFQBench =readBenchmark(base_demultiplexed+"/demultiplex_fq.benchmark")

################################################################################
#                           Compute Counts                                     #
################################################################################
if snakemake.config["gzip_input"] == "F":
    rawCounts = countFasta(base_sample+"/rawdata/fw.fastq", True);
else:
    rawCounts = countFastaGZ(base_sample+"/rawdata/fw.fastq.gz", True);
#rawCountsStr= '{0:g}'.format(float(rawCounts))
rawCountsStr= str(int(rawCounts))
#-peared
pearedCounts = 0
if snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T":
    pearedCounts = countFasta(base_sample+"/peared/seqs.assembled.fastq", True);
else:
    pearedCounts = countFasta(base_sample+"/peared/seqs.assembled.UNPAIRED.fastq", True);

#pearedCountsStr='{0:g}'.format(float(pearedCounts))
pearedCountsStr=str(int(pearedCounts))
prcPeared = "{:.2f}".format(float((pearedCounts/rawCounts)*100))
#-dumultiplex
if snakemake.config["demultiplexing"]["demultiplex"] != "F": #starting to test this  and snakemake.config["demultiplexing"]["bc_mismatch"]>0:
    #in the past we had two files fw and reverse nos everything is on one file
    #fwAssignedCounts = countFasta(base_sample_data+"/splitLibs/seqs.assigned.fna", False)
    #barcodes.fastq_corrected_toRC
    #rvAssignedCounts = countFasta(base_sample_data+"/splitLibsRC/seqs.assigned.fna", False)
    
    #prcFwAssigned = "{:.2f}".format(float((fwAssignedCounts/pearedCounts)*100))
    #prcRvAssigned = "{:.2f}".format(float((rvAssignedCounts/pearedCounts)*100))
    #totalAssigned = fwAssignedCounts + rvAssignedCounts
    #prcPearedAssigned = "{:.2f}".format(float((totalAssigned/pearedCounts)*100))
    #prcRawAssigned = "{:.2f}".format(float((totalAssigned/rawCounts)*100))
    #New implementation
    totalAssigned =  countFasta(base_split+"/seqs.assigned.fna", False)
    rvAssignedCounts = countTxt(base_sample+"/barcodes/barcodes.fastq_corrected_toRC")
    fwAssignedCounts = totalAssigned - rvAssignedCounts 
    prcFwAssigned = "{:.2f}".format(float((fwAssignedCounts/pearedCounts)*100))
    prcRvAssigned = "{:.2f}".format(float((rvAssignedCounts/pearedCounts)*100))
    prcPearedAssigned = "{:.2f}".format(float((totalAssigned/pearedCounts)*100))
    prcRawAssigned = "{:.2f}".format(float((totalAssigned/rawCounts)*100))

else: 
    totalAssigned = pearedCounts
    prcPearedAssigned = "{:.2f}".format(float((totalAssigned/pearedCounts)*100))
    prcRawAssigned = "{:.2f}".format(float((totalAssigned/rawCounts)*100))

#--cutadapt
cutSequences = False
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg":
    sequenceNoAdapters = countFasta(base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna", False)
    if (totalAssigned - sequenceNoAdapters) > 0:
        cutSequences = True
        prcCut = "{:.2f}".format(float((sequenceNoAdapters/totalAssigned)*100))
        prcCutRaw = "{:.2f}".format(float((sequenceNoAdapters/rawCounts)*100))

if removeChimeras:
    sequenceNoChimeras = countFasta(base_sample_data+"/seqs_fw_rev_filtered_nc.fasta", False)
    prcChim = "{:.2f}".format(float((sequenceNoChimeras/totalAssigned)*100))
    prcChimRaw = "{:.2f}".format(float((sequenceNoChimeras/rawCounts)*100))
    if cutSequences:
        prcChimCut = "{:.2f}".format(float((sequenceNoChimeras/sequenceNoAdapters)*100))
#out="{PROJECT}/runs/{run}/{sample}_data/"
trimmedCounts = countFasta(base_sample_data+"/seqs_fw_rev_filtered.fasta", False)
prcTrimmedSplit ="{:.2f}".format(float((trimmedCounts/totalAssigned)*100))
prcTrimmedRaw= "{:.2f}".format(float((trimmedCounts/rawCounts)*100))
if cutSequences:
    prcTrimmedCut="{:.2f}".format(float((trimmedCounts/sequenceNoAdapters)*100))
#if removeChimeras:
#    prcTrimmedChimera="{:.2f}".format(float((trimmedCounts/sequenceNoChimeras)*100))
try:
    samplesLib = subprocess.run( ["cat " + base_sample_data+"/seqs_fw_rev_filtered.dist.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    samplesLibInt = int(samplesLib.stdout.decode('utf-8').strip())
except Exception as e:
    totalReads = "Problem reading outputfile"
################################################################################
#                        Quality tool section                    #
################################################################################
qcStr = ""
if snakemake.config["QC"].lower() == "fastqc" or snakemake.config["QC"].lower() == "both":
    qcStr += ":red:`Tool:` [FastQC]_\n\n"
    qcStr += ":red:`Version:` "+ fqVersion +"\n\n"
    qcStr += "**Command:**\n\n"
    qcStr += ":commd:`fastqc "+ base_sample + "/rawdata/fw.fastq " + base_sample + "/rawdata/rv.fastq" + "--extract -o  "+ base_sample+"/qc/`\n\n"
    qcStr += "You can follow the links below, in order to see the complete FastQC report:\n\n"
    qcStr += "**Output files:**\n\n:green:`- FastQC for sample: "+snakemake.wildcards.sample+"_1:`  FQ1_ \n\n"
    qcStr += ".. _FQ1: ../../../samples/"+snakemake.wildcards.sample+"/qc/fw_fastqc.html \n\n"
    qcStr += "green:`- FastQC for sample: "+snakemake.wildcards.sample+"_2:`  FQ2_ \n\n"
    qcStr += ".. _FQ2: ../../../samples/"+snakemake.wildcards.sample+"/qc/rv_fastqc.html \n\n" 
    fqBench = readBenchmark(base_sample+"/qc/fq.benchmark")
    qcStr += fqBench + "\n\n"
if snakemake.config["QC"].lower() == "sequali" or snakemake.config["QC"].lower() == "both":
    qcStr += ":red:`Tool:` [Sequali]_\n\n"
    qcStr += ":red:`Version:` "+ sequaliVersion +"\n\n"
    qcStr += "**Command:**\n\n"
    qcStr += ":commd:`sequali --html  sequali.html --json sequali.json  --outdir  "+ base_sample+"/qc/ "+ base_sample + "/rawdata/fw.fastq " + base_sample + "/rawdata/rv.fastq`\n\n"
    qcStr += "You can follow the links below, in order to see the complete sequali report:\n\n"
    qcStr += "**Output files:**\n\n:green:`- Sequali report:`  SEQ1_ \n\n"
    qcStr += ".. _SEQ1: ../../../samples/"+snakemake.wildcards.sample+"/qc/sequali.html \n\n"
    qcStr += ":green:`- json report:` "+base_sample+"/qc/sequali.json \n\n"
    seqBench = readBenchmark(base_sample+"/qc/sequali.benchmark")
    qcStr += seqBench + "\n\n"    


################################################################################
#                         Generate sequence amounts chart                      #
################################################################################
numbers=[rawCounts,pearedCounts];
labels=["Raw", "Assembled"];
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    numbers.append(totalAssigned)
    labels.append("Demultiplexed")
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg":
    numbers.append(sequenceNoAdapters)
    labels.append("Cutadapt")
numbers.append(trimmedCounts)
labels.append("Length filtering")
if snakemake.config["chimera"]["search"] == "T" and removeChimeras:
    numbers.append(sequenceNoChimeras)
    labels.append("No Chimera")
createChart(numbers, tuple(labels),base_run+"/report_files/sequence_numbers."+snakemake.wildcards.sample+".png")
################################################################################
#                          Chimera check                                       #
################################################################################
variable_refs=""
if snakemake.config["chimera"]["search"] == "T" and snakemake.config["chimera"]["method"] == "usearch61":
    variable_refs+= ".. [usearch61] Edgar RC. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-2461.\n\n"
else: 
    variable_refs+= ".. [uchime] Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R (2011) UCHIME improves sensitivity and speed of chimera detection. Bioinformatics, 27 (16): 2194-2200. doi:10.1093/bioinformatics/btr381. \n\n"
quimeraStr = ""
if snakemake.config["chimera"]["search"] == "T":
    quimeraStr="Identify Chimera\n-------------------\n\n"
    quimeraStr+="Identify possible chimeric sequences (sequences generated due to the PCR amplification of multiple templates or parent sequences).\n\n"
    if snakemake.config["chimera"]["method"] == "usearch61":
        quimeraStr += ":red:`Tool:` [QIIME]_ - identify_chimeric_seqs.py\n\n"
        quimeraStr += ":red:`Version:` "+ icVersion +"\n\n"
        quimeraStr += ":red:`Method:` [usearch61]_ \n\n"
    else:
        quimeraStr += ":red:`Tool:` [Vsearch]_ - vsearch\n\n"
        quimeraStr += ":red:`Version:` "+ vsearchVersion +"\n\n"        
        quimeraStr += ":red:`Method:` "+ str(snakemake.config["chimera"]["method"]) +" - uses [uchime]_ \n\n"    
    quimeraStr += "**Command:**\n\n"
    if snakemake.config["chimera"]["method"] == "usearch61":
         quimeraStr+=":commd:`identify_chimeric_seqs.py -m "+ snakemake.config["chimera"]["method"]+" -i "+ base_sample_data+ "/seqs_fw_rev_accepted.fna "+str(snakemake.config["chimera"]["extra_params"])
         quimeraStr+=" -o "+ base_sample_data + "/chimera/` \n\n"
    else:
         quimeraStr+=":commd:`vsearch --"+ str(snakemake.config["chimera"]["method"])+" "+ base_sample_data + "/seqs_fw_rev_accepted.fna --threads "+ str(snakemake.config["chimera"]["threads"]) +" " +str(snakemake.config["chimera"]["extra_params"])   
         quimeraStr+=" --uchimeout "+ base_sample_data + "/chimera/chimeras.summary.txt` \n\n"
    quimeraStr+="**Output files:**\n\n"
    if snakemake.config["chimera"]["method"] == "usearch61":
        quimeraStr+=":green:`- File with the possible chimeric sequences:` "+base_sample_data + "/chimera/chimeras.txt\n\n"
    else:
        quimeraStr+=":green:`- File with the possible chimeric sequences:` "+base_sample_data + "/chimera/chimeras.summary.txt\n\n"
    identifyChimeraBench=readBenchmark(base_sample_data+"/chimera/chimera.benchmark")
    quimeraStr+=identifyChimeraBench
    quimeraStr+=chimera_log
    if removeChimeras:
        quimeraStr+=":red:`Reads after remove chimeric sequences:` "+ str(sequenceNoChimeras)+"\n\n"
        quimeraStr+=":red:`Percentage of reads vs raw reads:` "+ str(prcChimRaw) + "%\n\n"
        quimeraStr+=":red:`Percentage of reads vs demultiplexed reads:` "+ str(prcChim) + "%\n\n"
        if cutSequences:
            quimeraStr+=":red:`Percentage of reads vs cutadapt:` "+ str(prcChimRaw) + "%\n\n"




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
    fastQCPearStr += ":commd:`fastqc "+base_sample+"/peared/seqs.assembled.fastq --extract -o  "+ base_sample+"/peared/qc`\n\n"
    fastQCPearStr += "**Output files:**\n\n:green:`- FastQC report:` "+base_sample+"/peared/qc/seqs.assembled_fastqc.html FQ_Report_ \n\n"
    fastQCPearStr += ".. _FQ_Report: peared/qc/seqs.assembled_fastqc.html \n\n"
    fastQCPearStrBench =readBenchmark(base_sample+"/peared/qc/fq.benchmark")
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
    extractBCStr +=":commd:`extract_barcodes.py -f "+base_sample+"/peared/seqs.assembled.fastq -c "+str(snakemake.config["ext_bc"]["c"])+ " " + str(snakemake.config["ext_bc"]["bc_length"])+ " " + snakemake.config["ext_bc"]["extra_params"] + " -o "+base_barcodes+"/`\n\n"
    extractBCStr +="**Output files:**\n\n"
    extractBCStr +=":green:`- Fastq file with barcodes:` "+base_barcodes+"/barcodes.fastq\n\n"
    extractBCStr +=":green:`- Fastq file with the reads:` "+base_barcodes+"/reads.fastq\n\n"
    extractBCStr +=barBench
################################################################################
#                           CORRECT Barcodes                                   #
################################################################################
correctBCStr = ""
bcFile="barcodes.fastq"
if snakemake.config["demultiplexing"]["demultiplex"] != "F": # and snakemake.config["demultiplexing"]["bc_mismatch"]:
    correctBCStr = "Correct Barcodes\n--------------------\n"
    correctBCStr += "Try to correct the barcode from unassigned reads and place reads in correct orientetion.\n\n"
    correctBCStr += "Maximum number of mismatches **"  + str(snakemake.config["demultiplexing"]["bc_mismatch"]) + "**.\n\n"
    correctBCStr +=":red:`Tool:` Cascabel Java application\n\n"
    correctBCStr +="**Command:**\n\n"
    correctBCStr += ":commd:`java -jar Scripts/BarcodeCorrector.jar -b "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -fb "+base_barcodes+"/barcodes.fastq -fr "+ base_barcodes+"/reads.fastq  -m "  + str(snakemake.config["demultiplexing"]["bc_mismatch"]) + " -o  " + base_barcodes+"/barcodes.fastq_corrected -or  " + base_barcodes+"/reads.fastq_corrected -rc -x " + base_barcodes+"/sample_matrix.txt  >  " + base_barcodes+"/demux.log`\n\n"
    correctBCStr += "**Output files:**\n\n:green:`- Barcode corrected file:` "+base_barcodes+"/barcodes.fastq_corrected\n\n"
    correctBCStr += ":green:`- Reads corrected file:` "+base_barcodes+"/reads.fastq_corrected\n\n"
    correctBCStr += ":green:`- Error correction summary:` "+base_barcodes+"/demux.log\n\n"
    correctBarBench =readBenchmark(base_barcodes+"/barcodes_corrected.benchmark")
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
    splitStr+=":commd:`split_libraries_fastq.py -m "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -i "+base_barcodes+"/reads.fastq -o  "+base_split+" -b "+base_barcodes+"/"+bcFile+" -q "+str(snakemake.config["split"]["q"])+" -r "+str(snakemake.config["split"]["r"])+" --retain_unassigned_reads "+str(snakemake.config["split"]["extra_params"])+" --barcode_type "+str(snakemake.config["split"]["barcode_type"])+"`\n\n"
    splitStr+=splitLibsBench

    splitStr+="Retain assigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`cat "+base_split+"/seqs.fna | grep -P -A1 \"(?!>Unass)^>\" | sed '/^--$/d' > "+base_split+"/seqs.assigned.fna`\n\n"

    splitStr+="Create file with only unassigned reads\n"
    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    splitStr+="**Command:**\n\n"
    splitStr+=":commd:`cat "+base_split+"/seqs.fna | grep \"^>Unassigned\" |  sed 's/>Unassigned_[0-9]* /@/g' | sed 's/ .*//' | grep -F -w -A3  -f - "+base_sample+"/peared/seqs.assembled.fastq |  sed '/^--$/d' >"+base_split+"/unassigned.fastq`\n\n"

#    splitStr+="Reverse complement unassigned reads\n"
#    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
#    splitStr+=":red:`Tool:` [Vsearch]_\n\n"
#    splitStr+=":red:`version:`  "+vsearchVersion+"\n\n"
#    splitStr+="**Command:**\n\n"
#    splitStr+=":commd:`vsearch --fastx_revcomp "+base_sample_data+"/splitLibs/unassigned.fastq  --fastqout "+base_sample_data+"/splitLibs/unassigned.reversed.fastq`\n\n"


#    splitStr+="Barcode extraction for reverse complemented, unassigned reads\n"
#    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
#    splitStr +=":red:`Tool:` [QIIME]_ - extract_barcodes.py\n\n"
#    splitStr +=":red:`Version:` "+ebVersion+"\n\n"
#    splitStr+="**Command:**\n\n"
#    splitStr +=":commd:`extract_barcodes.py -f "+base_sample_data+"/splitLibs/unassigned.reversed.fastq -c "+str(snakemake.config["ext_bc"]["c"])+" "+str(snakemake.config["ext_bc"]["bc_length"])+" "+snakemake.config["ext_bc"]["extra_params"]+" -o "+base_sample_data+"/barcodes_unassigned/`\n\n"

#    if snakemake.config["bc_mismatch"]:
#        splitStr += "Correct reverse complemented barcodes \n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
#        splitStr += "Maximum number of mismatches **"  + str(snakemake.config["bc_mismatch"]) + "**.\n\n"
#        splitStr +=":red:`Tool:` Cascabel Java application\n\n"
#        splitStr +="**Command:**\n\n"
#        splitStr += ":commd:`java -cp Scripts/BarcodeCorrector/build/classes/  barcodecorrector.BarcodeCorrector -b "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -f "+base_sample_data+"/barcodes_unassigned/barcodes.fastq_corrected -m "  + str(snakemake.config["bc_mismatch"]) + "`\n\n"
#        splitStr += "**Output file:**\n\n:green:`- Barcode corrected file:` "+base_barcodes+"/barcodes.fastq_corrected\n\n"
#        splitStrBench =readBenchmark(base_sample_data+"/barcodes_unassigned/barcodes_corrected.benchmark")
#        splitStr += splitStrBench+"\n\n"

#    splitStr +="Split reverse complemented reads\n"
#    splitStr+="~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
#    splitStr +=":red:`Tool:` [QIIME]_ - extract_barcodes.py\n\n"
#    splitStr +=":red:`Version:` "+ebVersion+"\n\n"
#    splitStr+="**Command:**\n\n"
#    splitStr +=":commd:`split_libraries_fastq.py -m "+snakemake.wildcards.PROJECT+"/metadata/sampleList_mergedBarcodes_"+snakemake.wildcards.sample+".txt -i "+base_sample_data+"/barcodes_unassigned/reads.fastq -o "+base_sample_data+"/splitLibsRC -b "+base_run+"/"+str(snakemake.wildcards.sample)+"_data/barcodes_unassigned/"+bcFile+" -q "+str(snakemake.config["split"]["q"])+" -r "+str(snakemake.config["split"]["r"])+" "+str(snakemake.config["split"]["extra_params"])+" --barcode_type "+str(snakemake.config["split"]["barcode_type"])+"`\n\n"
#    splitStr +=splitLibsBench+"\n\n"

    splitStr +="**Output files:**\n\n"
#   # splitStr +=":green:`- FW reads fasta file with new header:` "+base_sample_data+"/splitLibs/seqs.assigned.fna\n\n"
    splitStr +=":green:`- Text histogram with the length of the fw reads:` "+base_split+"/histograms.txt\n\n"
    splitStr +=":green:`- Log file for the fw reads:` "+base_split+"/split_library_log.txt\n\n"
#   # splitStr +=":green:`- RV reads fasta file with new header:` "+base_sample_data+"/splitLibsRC/seqs.assigned.fna\n\n"
#    splitStr +=":green:`- Text histogram with the length of the rv reads:` "+base_sample_data+"/splitLibsRC/histograms.txt\n\n"
#    splitStr +=":green:`- Log file for the rv reads:` "+base_sample_data+"/splitLibsRC/split_library_log.txt\n\n"
#    splitStr +=":green:`- Fasta file with unassigned reads:` "+base_sample_data+"/splitLibsRC/seqs.unassigned.fna\n\n"
    splitStr +=":red:`Number of reads assigned on FW:` "+str(fwAssignedCounts)+" = "+str(prcFwAssigned)+"% of the peared reads\n\n"
    splitStr +=":red:`Number of reads assigned on RVC:` "+str(rvAssignedCounts)+" = "+str(prcRvAssigned)+"% of the peared reads\n\n"

################################################################################
#                           Single FastQ creation                              #
################################################################################
demultiplexFQ = ""
if snakemake.config["demultiplexing"]["demultiplex"] == "T" and snakemake.config["demultiplexing"]["create_fastq_files"] == "T":
    demultiplexFQ = "Generate single sample fastq files\n------------------------------------------\n\n" # title
    demultiplexFQ += "Create single fastq files per samples (based on the raw data without applying any filtering).\n\n"
    demultiplexFQ +=":red:`Tool:` Cascabel Java program\n\n"
    demultiplexFQ += "**Command:**\n\n"
    demultiplexFQ += ":commd:`"+snakemake.config["java"]["command"] + " -cp Scripts DemultiplexQiime --txt -a rv -b "+ str(snakemake.config["demultiplexing"]["bc_mismatch"]) + " -d "+ base_split+"/seqs.assigned.ori.txt -o "+ base_demultiplexed+"/ "
    ext=".gz"
    if snakemake.config["gzip_input"].casefold() == "f":
        ext=""
    demultiplexFQ += "-r1 "+base_sample+"/rawdata/fw.fastq"+ext+" -r2 "+base_sample+"/rawdata/fw.fastq"+ext+"`\n\n"
    if snakemake.config["demultiplexing"]["remove_bc"]:
        demultiplexFQ +=":red:`Barcodes removed:` "+ str(snakemake.config["demultiplexing"]["remove_bc"]) + " first bases\n\n"
#Now only for ASV workflow
   # if snakemake.config["primers"]["remove"].lower() == "cfg":
   #     demultiplexFQ +=":red:`Primers removed:` **FW** " + snakemake.config["primers"]["fw_primer"] + " **RV** " +snakemake.config["primers"]["rv_primer"]+"\n\n"
   # elif snakemake.config["primers"]["remove"].lower() == "metadata":
   #     demultiplexFQ +=":red:`Removed primers` were obtained from the metadata file.\n\n" 
    demultiplexFQ += "**The demultiplexed fastq files are located at:**\n\n:green:`- Demultiplexed directory:` "+base_demultiplexed+"/\n\n"
    demultiplexFQ += ":green:`- Summary file:` "+base_demultiplexed+"/summary.pcr.txt\n\n"
    demultiplexFQ += demultiplexFQBench
   # also this only for the ASV workflow
   # if (snakemake.config["primers"]["remove"].lower() == "cfg" or snakemake.config["primers"]["remove"].lower() == "metadata"):
   #     demultiplexFQ += "**Remove primers:**\n\nFollowing, primers were removed from the fastq files\n\n"
   #     demultiplexFQ +=":red:`Tool:` [Cutadapt]_\n\n"
   #     demultiplexFQ += ":red:`Version:` "+cutVersion+"\n\n"
   #     demultiplexFQ += "**Command:**\n\n"
   #     if snakemake.config["primers"]["remove"].lower() == "cfg":
   #         if snakemake.config["LIBRARY_LAYOUT"].casefold()=="pe":
   #             demultiplexFQ += ":commd:`cutadapt -g "+ snakemake.config["primers"]["fw_primer"]  + " -G " + snakemake.config["primers"]["rv_primer"]  + " " +snakemake.config["primers"]["extra_params"]+" -O "+ snakemake.config["primers"]["min_overlap"]  +" -m " +snakemake.config["primers"]["min_length"]+ " -o "+base_demultiplexed+"/primer_removed/SAMPLE_1.fastq.gz -p "+base_demultiplexed+"/primer_removed/SAMPLE_2.fastq.gz "+ base_demultiplexed+"/SAMPLE_1.fq.gz  "+ base_demultiplexed+"/SAMPLE_2.fq.gz  >> "+base_demultiplexed+"/primer_removed/"+snakemake.wildcards.sample+".cutadapt.log`\n\n"
   #         else:
   #             demultiplexFQ += ":commd:`cutadapt -g "+ snakemake.config["primers"]["fw_primer"]  + " " +snakemake.config["primers"]["extra_params"]+" -O "+ snakemake.config["primers"]["min_overlap"]  +" -m " +snakemake.config["primers"]["min_length"]+ " -o "+base_demultiplexed+"/primer_removed/SAMPLE_1.fastq.gz "+ base_demultiplexed+"/SAMPLE_1.fq.gz  >> "+base_demultiplexed+"/primer_removed/"+snakemake.wildcards.sample+".cutadapt.log`\n\n"
   #         demultiplexFQ += "The above command ran once for each single sample fastq file(s) using the mentioned primers\n\n"
   #     else: #is from metadata
   #         if snakemake.config["LIBRARY_LAYOUT"].casefold()=="pe":
   #             demultiplexFQ += ":commd:`cutadapt -g sample_FW_primer  -G sample_RV_primer " +snakemake.config["primers"]["extra_params"]+" -O "+ snakemake.config["primers"]["min_overlap"]  +" -m " +snakemake.config["primers"]["min_length"]+ " -o "+base_demultiplexed+"/primer_removed/SAMPLE_1.fastq.gz -p "+base_demultiplexed+"/primer_removed/SAMPLE_2.fastq.gz "+ base_demultiplexed+"/SAMPLE_1.fq.gz "+ base_demultiplexed+"/SAMPLE_2.fq.gz  >> "+base_demultiplexed+"/primer_removed/"+snakemake.wildcards.sample+".cutadapt.log`\n\n"
   #         elif snakemake.config["LIBRARY_LAYOUT"].casefold()=="se":
   #             demultiplexFQ += ":commd:`cutadapt -g sample_FW_primer "+ " " +snakemake.config["primers"]["extra_params"]+" -O "+ snakemake.config["primers"]["min_overlap"]  +" -m " +snakemake.config["primers"]["min_length"]+ " -o "+base_demultiplexed+"/primer_removed/SAMPLE_1.fastq.gz "+ base_demultiplexed+"/SAMPLE_1.fq.gz  >> "+base_demultiplexed+"/primer_removed/"+snakemake.wildcards.sample+".cutadapt.log`\n\n"
   #         demultiplexFQ += "The above command ran once for each single sample fastq file(s) and primers were obtained from the mapping file accordingly to its sample\n\n"    
   #     demultiplexFQ += ":green:`- Reads without primers:` "+base_demultiplexed+"/primer_removed\n\n"
   #     if "--discard-untrimmed" in snakemake.config["primers"]["extra_params"]:
   #         demultiplexFQ += ":green:`- Discarded reads (no primer):` "+base_demultiplexed+"/reads_discarded_primer\n\n"
   #     else:
   #         demultiplexFQ += ":red:`- Given the options, reads without primers where not removed!`\n\n"
   #     demultiplexFQ += ":green:`- Primer removal results by sample:` primers_removal_\n\n"
   #     demultiplexFQ +=" .. _primers_removal: report_files/cutadapt."+snakemake.wildcards.sample+".fastq_summary.tsv\n\n"

################################################################################
#                           Combine FW and Reverse reads                       #
################################################################################

combineFR = ""
#if snakemake.config["demultiplexing"]["demultiplex"] != "F":
#    combineFR = "Combine reads\n---------------------------------\n\n" # title
#    combineFR += "Concatenate forward and reverse reads.\n\n"
#    combineFR += "**Command:**\n\n"
#    combineFR += ":commd:`cat "+base_sample_data+"/splitLibs/seqs.assigned.fna "+base_sample_data+"/splitLibsRC/seqs.assigned.fna > "+base_sample_data+"/seqs_fw_rev_accepted.fna`\n\n"
#    combineFR +="**Output files:**\n\n"
#    combineFR +=":green:`- Fasta file with combined reads:` "+base_sample_data+"/seqs_fw_rev_accepted.fna\n\n"
#    combineFR +=":red:`- Total number of acepted reads:` " +str(totalAssigned)+ " = "+ str(prcPearedAssigned)+ "% of the peared reads or "+str(prcRawAssigned)+"% of the raw reads.\n\n"
#    combineFR += combineBench
    
################################################################################
#                          Cut adapters                                        #
################################################################################
cutAdaptStr = ""
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg":
    cutAdaptStr = "Remove sequence primers\n------------------------\n\n" # title
    cutAdaptStr +="Remove the adapters / primers from the reads.\n\n"
    cutAdaptStr +=":red:`Tool:` [Cutadapt]_\n\n"
    cutAdaptStr += ":red:`Version:` "+cutVersion+"\n\n"
    cutAdaptStr += "**Command:**\n\n"
    primer_lines=0
    if snakemake.config["primers"]["remove"].lower() == "cfg":
        #cutAdaptStr += ":commd:`cutadapt "+ str(snakemake.config["cutadapt"]["adapters"])+" " + str(snakemake.config["cutadapt"]["extra_params"]) + " -o " + base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna\n\n"
        #cutAdaptStr +=  base_sample_data+"/seqs_fw_rev_accepted.fna > " +  base_sample_data+"/seqs_fw_rev_accepted_no_adapters.log`\n\n"
        cutAdaptStr += ":commd:`cutadapt -g "+ str(snakemake.config["primers"]["fw_primer"])+"..."+str(snakemake.config["primers"]["rv_primer"])+" "+ str(snakemake.config["primers"]["extra_params"]) + " -o " + base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna "+base_sample_data+"/seqs_fw_rev_accepted.fna > " +  base_sample_data+"/seqs_fw_rev_accepted_no_adapters.log`\n\n"
    elif snakemake.config["primers"]["remove"].lower() == "metadata":
        primers=""
        try:
            #with open(base_sample_data+"/primers.txt") as pfile:
            with open(base_run+"/report_files/primers."+snakemake.wildcards.sample+".txt") as pfile:
                primers=pfile.read()
                #primer_lines=len(pfile.readlines())
                primer_lines=len(primers.split("\n"))
                if primer_lines > 1:
                    if snakemake.config["LIBRARY_LAYOUT"].casefold()=="pe":
                        primers="-g sample_FW_primer...sampleRV_primer"
                    else:
                        primers="-g sample_FW_primer"

        except FileNotFoundError:
            primers="-ERROR reading primer file-"
        #cutAdaptStr += ":commd:`cutadapt "+primers +" " + str(snakemake.config["cutadapt"]["extra_params"]) + " -o " + base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna\n\n"
        #cutAdaptStr += base_sample_data+"/seqs_fw_rev_accepted.fna > " +  base_sample_data+"/seqs_fw_rev_accepted_no_adapters.log`\n\n"
        cutAdaptStr += ":commd:`cutadapt "+primers +" " + str(snakemake.config["primers"]["extra_params"]) + " -o " + base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna "+  base_sample_data+"/seqs_fw_rev_accepted.fna > " +  base_sample_data+"/seqs_fw_rev_accepted_no_adapters.log`\n\n"
        
#cutAdaptStr += "*PRIMERS: primer sequences were obtained from the metadata file\n\n"
    if primer_lines > 1:
        cutAdaptStr += ":green:`- Primers used by sample:` primers_sample_\n\n"
        cutAdaptStr +=  ".. _primers_sample: report_files/primers."+snakemake.wildcards.sample+".txt\n\n"
    cutAdaptStr += "**Output files:**\n\n:green:`- Reads without adapters:` "+base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna\n\n"
    if cutSequences:
        cutAdaptStr += ":red:`Total number of reads after cutadapt:` "+ str(sequenceNoAdapters) + " = " + str(prcCut) + "% of the assigned reads or "+ str(prcCutRaw)+"% of the total reads\n\n"
    #cutAdaptStr+=":\n\n"
    cutAdaptStr+=":green:`- Primer removal results by sample:` primers_OTU_\n\n"
    cutAdaptStr+=" .. _primers_OTU: report_files/cutadapt."+snakemake.wildcards.sample+".summary.tsv\n\n"

    cutAdaptBench =readBenchmark(base_sample_data+"/cutadapt.benchmark")
    cutAdaptStr += cutAdaptBench+"\n\n"
################################################################################
#                          Counts for too long too shorts                      #
################################################################################
#trimmedStr =  ":red:`Total number of reads after trimming:` "+str(trimmedCounts)+ "="+ str(prcTrimmedSplit)+"% of the demultiplexed reads or " + str(prcTrimmedRaw) + "% of the raw reads\n\n"
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
shorts = str(snakemake.config["rm_reads"]["shorts"])
longs = str(snakemake.config["rm_reads"]["longs"])
with open(base_sample_data+"/filter.log") as trimlog:
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
data.append(base_sample+"/rawdata/\*.fq")
data.append(str(rawCounts))
data.append("{:.2f}".format(float((rawCounts/rawCounts)*100))+"%")
fileData.append(data)
data=[]
#pear
data.append("Assembled reads")
data.append(base_sample+"/peared/seqs.assembled.fastq")
data.append(str(pearedCounts))
data.append("{:.2f}".format(float((pearedCounts/rawCounts)*100))+"%")
fileData.append(data)
data=[]
#splitted
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    data.append("Demultiplexed reads")
    data.append(base_sample_data+"/seqs_fw_rev_accepted.fna")
    data.append(str(totalAssigned))
    data.append("{:.2f}".format(float((totalAssigned/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
#adapters
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg":
    data.append("Adapter removed")
    data.append(base_sample_data+"/seqs_fw_rev_accepted_no_adapters.fna")
    data.append(str(sequenceNoAdapters))
    data.append("{:.2f}".format(float((sequenceNoAdapters/rawCounts)*100))+"%")
    fileData.append(data)
    data=[]
#length filtered
data.append("Length filtered")
data.append(base_sample_data+"/seqs_fw_rev_filtered.fasta")
data.append(str(trimmedCounts))
data.append("{:.2f}".format(float((trimmedCounts/rawCounts)*100))+"%")
fileData.append(data)
data=[]
#chimera
if snakemake.config["chimera"]["search"] == "T" and removeChimeras:
    data.append("Non chimeric reads")
    data.append(base_sample_data+"/seqs_fw_rev_filtered_nc.fasta")
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
    dist_table = readSampleDist(base_sample_data+"/seqs_fw_rev_filtered.dist.txt",trimmedCounts,samplesLibInt)
    sampleDistChart = "Sample distribution\n--------------------------------------\n\n" # title
    sampleDistChart += dist_table + "\n\n"
    sampleDistChart += ".. image:: report_files/seqs_fw_rev_filtered."+snakemake.wildcards.sample+".dist.png\n\n"
    sampleDistChart +="The previous chart shows the number of clean reads per sample. The bars are sorted from left to right, according to the metadata input file.\n\n"
    sampleDistChart +="**To see more details about the number of reads per sample in this library, please refer to the file:** "+base_sample_data+"/seqs_fw_rev_filtered.dist.txt\n\n"


################################################################################
#                       User description section                               #
################################################################################
desc = snakemake.config["description"]
txtDescription = ""
if len(desc) > 0:
    txtDescription = "\n**User description:** "+desc+"\n"


################################################################################
#                       controls warning section                               #
################################################################################
"""
We want to include a small section to warn the user about the use of controls. This could be
the case if they are demultiplexing a complete library. 
"""
ctrlWarning =""
if snakemake.config["demultiplexing"]["demultiplex"] == "T":
    ctrlWarning="\n:warn:`Note: Library demultiplexing has been carried out, if you have controls among your samples, please be aware that Cascabel won't perform any special operation with them. They are treated as any other sample within this workflow. Please make sure to analyze your controls with other tools, and correct your sample counts for potential contamination.`\n"
################################################################################
#                                Report                                        #
################################################################################
taxoTool = str(snakemake.config["assignTaxonomy"]["tool"].lower())

report("""
Amplicon Analysis Report for Library: {snakemake.wildcards.sample}
=====================================================================
    .. role:: commd
    .. role:: red
    .. role:: green
    .. role:: warn

**CASCABEL** is designed to run amplicon sequence analysis across single or multiple read libraries.

The objective of this pipeline is to create different output files which allow the user to explore data in a simple and meaningful way, as well as facilitate downstream analysis, based on the generated output files.

Another aim of **CASCABEL** is also to encourage the documentation process, by creating this report in order to assure data analysis reproducibility.

:red:`Cascabel version:` {version}

{txtDescription}

{ctrlWarning}

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

{qcStr}


Read pairing
----------------
Align paired end reads and merge them into one single sequence in case they overlap.

:red:`Tool:` [PEAR]_

:red:`version:` {pearversion}

**Command:**

:commd:`pear -f {base_sample}/rawdata/fw.fastq -r {base_sample}/rawdata/rv.fastq -t {snakemake.config[pear][t]} -v {snakemake.config[pear][v]} -j {snakemake.config[pear][j]} -p {snakemake.config[pear][p]} -o {base_sample_data}/peared/seqs > {base_sample_data}/peared/seqs.assembled.fastq`

**Output files:**

:green:`- Merged reads:` {base_sample_data}/peared/seqs.assembled.fastq

:green:`- Log file:` {base_sample_data}/peared/pear.log

:red:`Number of peared reads:` {pearedCountsStr} =  {prcPeared}%

{pearBench}

{fastQCPearStr}

{extractBCStr}

{correctBCStr}

{splitStr}

{demultiplexFQ}

{combineFR}

{cutAdaptStr}


Remove too long and too short reads
------------------------------------
Remove very short and long reads, with lengths more than some standard deviation below or above the mean to be short or long respectively

:green:`- Minimun length expected (shorts):` {shorts}

:green:`- Maximun length expected (longs):` {longs}

**Command:**

:commd:`awk '!/^>/ {{ next }} {{ getline seq }} length(seq) > shorts  && length(seq) < longs {{ print $0 \"\\n\" seq }}'  {base_sample_data}/seqs_fw_rev_accepted.fna  >  {base_sample_data}/seqs_fw_rev_filtered.fasta`

**Sequence distribution before remove reads**

.. image:: report_files/seqs_dist_hist.{snakemake.wildcards.sample}.png
    :height: 400px
    :width: 400px
    :align: center


**Output file:**

:green:`- Fasta file with correct sequence length:` {base_sample_data}/seqs_fw_rev_filtered.fasta

{trimmedStr}

{rmShorLongBench}


{quimeraStr}


{sampleDistChart}


Final counts
-------------

{countTxt}

.. image:: report_files/sequence_numbers.{snakemake.wildcards.sample}.png

OTU report
---------------------------

Cascabel report on downstream analyses in combination with multiple libraries (if supplied), can be found at the following link: otu_report_ ({base_run}/otu_report_{taxoTool}.html)

    .. _otu_report: otu_report_{taxoTool}.html

References
------------------

.. [FastQC] FastQC v0.11.3. Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data

.. [PEAR] PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Zhang et al (2014) Bioinformatics 30(5): 614-620 | doi:10.1093/bioinformatics/btt593

.. [QIIME] QIIME. Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Gonzalez Pena A, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. 2010. QIIME allows analysis of high-throughput community sequencing data. Nature Methods 7(5): 335-336.

.. [Cutadapt] Cutadapt v1.15 .Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

.. [Vsearch] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584


{variable_refs}


""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
