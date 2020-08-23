import subprocess
import functools
from snakemake.utils import report
from benchmark_utils import readBenchmark
from benchmark_utils import countTxt
from seqsChart import createChart
from seqsChart import createChartPrc
from benchmark_utils import countFasta
from benchmark_utils import make_table

################
#Function to retrive the sample names and put in the report title
#@param file with the sample list, it is created during combine_filtered_samples
#snakemake.wildcards.project + "/runs/" + snakemake.wildcards.run + "/samples.log"
#@return the title with the samples
def getSampleList(sampleFile):
    with open(sampleFile) as sfile:
        samps ="Amplicon Analysis Report for Libraries: "
        for l in sfile:
            samps+= l
        samps+="\n"
        for i in range(0,len(samps)):
            samps+="="
    return samps;

#########################
#This function reads the file cat_samples.log which have the executed command to
#combine all the libraries after cleaning and demultiplexing and before taxonomy
#assignation
#@param catLogFile file with the command
#snakemake.wildcards.project + "/runs/" + snakemake.wildcards.run + "/cat_samples.log"
#@return the string ready to be concatenated into the report.
def getCombinedSamplesList(catLogFile):
    with open(catLogFile) as sfile:
        command =":commd:`"
        i=0
        for l in sfile:
            if i == 0:
                command+= l + "`\n\n"
            i+=1
    return command;


#title = getSampleList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/samples.log")
#catCommand =  getCombinedSamplesList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/cat_samples.log")
title = "Amplicon Analysis Report\n===========================\n\n"
################################################################################
#                         Benchmark Section                                    #
# This section is to generate a pre-formatted text with the benchmark info for #
# All the executed rules.                                                      #
################################################################################
#combineBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/combine_seqs_fw_rev.benchmark")
dada2Benchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2.benchmark")
asvFilterBenchmark =  readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/filter.benchmark")

#pikRepBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/pick_reps.benchmark")
#assignTaxaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/assign_taxa.benchmark")
otuTableBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/dada2.table.benchmark")
convertOtuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/dada2.biom.benchmark")
#convertOtuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable.txt.benchmark")
summTaxaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/summary/summarize_taxa.benchmark")
asvNoSingletonsBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/asvTable_nosingletons.bio.benchmark")
filterASVTableBenchmark =  readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/asvTable_nosingletons.txt.benchmark")
filterBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/representative_seq_set_noSingletons.benchmark")
deRepBenchmark=""
#if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
#    deRepBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/derep.benchmark")
if snakemake.config["alignRep"]["align"] == "T":
    #align_seqs.py -m {config[alignRep][m]} -i {input} -o {params.outdir} {config[alignRep][extra_params]}
    alignBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/align_rep_seqs.benchmark")
    #"filter_alignment.py -i {input} -o {params.outdir} {config[filterAlignment][extra_params]}"
    alignFilteredBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/filtered/align_rep_seqs.benchmark")
    #"make_phylogeny.py -i {input} -o {output} {config[makeTree][extra_params]}"
    makePhyloBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.benchmark")
kronaBenchmark=""
if snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kronaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/krona_report.benchmark")

#dada2FilterBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/filter.benchmark")
#dada2Benchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2.benchmark")
#dada2BiomBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2.biom.benchmark")

################################################################################
#                         TOOLS VERSION SECTION                          #
################################################################################
#clusterOtuV = subprocess.run([snakemake.config["qiime"]["path"]+'pick_otus.py', '--version'], stdout=subprocess.PIPE)
#clusterOtuVersion = "**" + clusterOtuV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

#pickRepV = subprocess.run([snakemake.config["qiime"]["path"]+'pick_rep_set.py', '--version'], stdout=subprocess.PIPE)
#pickRepVersion = "**" + pickRepV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

#assignTaxaV = subprocess.run([snakemake.config["qiime"]["path"]+'parallel_assign_taxonomy_'+snakemake.config["assignTaxonomy"]["qiime"]["method"]+'.py', '--version'], stdout=subprocess.PIPE)
#assignTaxaVersion = "**" + assignTaxaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

#makeOTUV = subprocess.run([snakemake.config["qiime"]["path"]+'make_otu_table.py', '--version'], stdout=subprocess.PIPE)
#makeOTUVersion = "**" + makeOTUV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

convertBiomV = subprocess.run([snakemake.config["biom"]["command"], '--version'], stdout=subprocess.PIPE)
convertBiomVersion = "**" + convertBiomV.stdout.decode('utf-8').strip() + "**"

dada2V = subprocess.run([snakemake.config["Rscript"]["command"],'Scripts/dada2Version.R'], stdout=subprocess.PIPE)
dada2Version = "**" + dada2V.stdout.decode('utf-8').strip() + "**"


summTaxaSV = subprocess.run([snakemake.config["qiime"]["path"]+'summarize_taxa.py', '--version'], stdout=subprocess.PIPE)
summTaxaVersion = "**" + summTaxaSV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

filterOTUNoSV = subprocess.run([snakemake.config["qiime"]["path"]+'filter_otus_from_otu_table.py', '--version'], stdout=subprocess.PIPE)
filterOTUNoSVersion = "**" + filterOTUNoSV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

filterFastaV = subprocess.run([snakemake.config["qiime"]["path"]+'filter_fasta.py', '--version'], stdout=subprocess.PIPE)
filterFastaVersion = "**" + filterFastaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

rscriptV = subprocess.run([snakemake.config["Rscript"]["command"], '--version'], stdout=subprocess.PIPE)
rscriptVersion = "**" + filterFastaV.stdout.decode('utf-8').strip() + "**"


#blastnV = subprocess.run([snakemake.config["assignTaxonomy"]["blast"]["command"], '-version'], stdout=subprocess.PIPE)
#blastnVersion = "**" + blastnV.stdout.decode('utf-8').split('\n', 1)[0].replace('blastn:','').strip() + "**"

#vsearchV2 = subprocess.run([snakemake.config["assignTaxonomy"]["vsearch"]["command"], '--version'], stdout=subprocess.PIPE)
#vsearchVersion_tax = "**" + vsearchV2.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"

#if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
#    vsearchV = subprocess.run([snakemake.config["derep"]["vsearch_cmd"], '--version'], stdout=subprocess.PIPE)
#    vsearchVersion = "**" + vsearchV.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"

if snakemake.config["alignRep"]["align"] == "T":
    alignFastaVersion="TBD"
    try:
        alignFastaV = subprocess.run([snakemake.config["qiime"]["path"]+'align_seqs.py', '--version'], stdout=subprocess.PIPE)
        if "Version" in alignFastaVersion:
            alignFastaVersion = "**" + alignFastaV.stdout.decode('utf-8').replace('Version: ','').strip() + "**"
    except Exception as e:
        alignFastaVersion = "**Problem retriving the version**"

    filterAlignmentV = subprocess.run([snakemake.config["qiime"]["path"]+'filter_alignment.py', '--version'], stdout=subprocess.PIPE)
    filterAlignmentVersion = "**" + filterAlignmentV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

    makePhyloV = subprocess.run([snakemake.config["qiime"]["path"]+'make_phylogeny.py', '--version'], stdout=subprocess.PIPE)
    makePhyloVersion = "**" + makePhyloV.stdout.decode('utf-8').replace('Version:','').strip() + "**"


################################################################################
#                        Compute counts section                                #
################################################################################
totalReads = "TBD"
intTotalReads = 1;
try:
     treads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/filter_summary.out | awk 'NR>1{sum=sum+$2} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intTotalReads = int(treads.stdout.decode('utf-8').strip())
     totalReads = "**" + str(intTotalReads) + "**"
except Exception as e:
     totalReads = "Problem reading outputfile"

filteredReads = "TBD"
intFilteredReads = 1;
prcFiltered=0.0
try:
     freads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/filter_summary.out | awk 'NR>1{sum=sum+$3} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intFilteredReads = int(freads.stdout.decode('utf-8').strip())
     filteredReads = "**" + str(intFilteredReads) + "**"
     prcFiltered = float(intFilteredReads/intTotalReads)*100
     prcFilteredStr = "**" + "{:.2f}".format(prcFiltered) + "%**"
except Exception as e:
     filteredReads = "Problem reading outputfile"

denoisedFWReads = "TBD"
intDenoisedFWReads = 1;
prcDenoisedFW=0
try:
     dfwreads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/stats_dada2.txt | awk 'NR>1{sum=sum+$2} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intDenoisedFWReads = int(dfwreads.stdout.decode('utf-8').strip())
     denoisedFWReads = "**" + str(intDenoisedFWReads) + "**"
     prcDenoisedFW = float(intDenoisedFWReads/intTotalReads)*100
     prcDenoisedFWStr = "**" + "{:.2f}".format(prcDenoisedFW) + "%**"
     prcDenoisedFWvsFiltered = (intDenoisedFWReads/intFilteredReads)*100
     prcDenoisedFWStrvsFiltered = "**" + "{:.2f}".format(prcDenoisedFWvsFiltered) + "%**"
except Exception as e:
     denoisedFWReads = "Problem reading outputfile"

denoisedRVReads = "TBD"
intDenoisedRVReads = 1;
prcDenoisedRV=0.0
try:
     drvreads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/stats_dada2.txt | awk 'NR>1{sum=sum+$3} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intDenoisedRVReads = int(drvreads.stdout.decode('utf-8').strip())
     denoisedRVReads = "**" + str(intDenoisedRVReads) + "**"
     prcDenoisedRV = float(intDenoisedRVReads/intTotalReads)*100
     prcDenoisedRVStr = "**" + "{:.2f}".format(prcDenoisedRV) + "%**"
     prcDenoisedRVvsFiltered = (intDenoisedRVReads/intFilteredReads)*100
     prcDenoisedRVStrvsFiltered = "**" + "{:.2f}".format(prcDenoisedRVvsFiltered) + "%**"
except Exception as e:
     denoisedRVReads = "Problem reading outputfile"

mergedReads = "TBD"
intmergedReads = 1;
prcmerged=0.0
try:
     mreads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/stats_dada2.txt | awk 'NR>1{sum=sum+$4} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intmergedReads = int(mreads.stdout.decode('utf-8').strip())
     mergedReads = "**" + str(intmergedReads) + "**"
     prcmerged = float(intmergedReads/intTotalReads)*100
     prcmergedStr = "**" + "{:.2f}".format(prcmerged) + "%**"
     prcmergedvsVariant = (intmergedReads/((intDenoisedFWReads+intDenoisedFWReads)/2))*100
     prcmergedStrvsVariant = "**" + "{:.2f}".format(prcmergedvsFiltered) + "%**"
except Exception as e:
     mergedReads = "Problem reading outputfile"

lengthFReads = "TBD"
intlengthFReads = 1;
prclengthF=0.0
try:
     lreads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/stats_dada2.txt | awk 'NR>1{sum=sum+$5} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
     intlengthFReads = int(lreads.stdout.decode('utf-8').strip())
     lengthFReads = "**" + str(intlengthFReads) + "**"
     prclengthF = float(intlengthFReads/intTotalReads)*100
     prclengthFStr = "**" + "{:.2f}".format(prclengthF) + "%**"
     prclengthFvsMerged = (intlengthFReads/intmergedReads)*100
     prclengthFStrvsMerged = "**" + "{:.2f}".format(prclengthFvsMerged) + "%**"
except Exception as e:
     lengthFReads = "Problem reading outputfile"

chimeraReads = "TBD"
intchimeraReads = 1;
prcchimera=0.0
if snakemake.config["dada2_asv"]["chimeras"] == "T":
    try:
         chreads = subprocess.run( ["cat " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/stats_dada2.txt | awk 'NR>1{sum=sum+$6} END{print sum}'"], stdout=subprocess.PIPE, shell=True)    
         intchimeraReads = int(chreads.stdout.decode('utf-8').strip())
         chimeraReads = "**" + str(intchimeraReads) + "**"
         prcchimera = float(intchimeraReads/intTotalReads)*100
         prcchimeraStr = "**" + "{:.2f}".format(prcchimera) + "%**"
         prcchimeravsLength = (intchimeraReads/intlengthFReads)*100
         prcchimeraStrvsLength = "**" + "{:.2f}".format(prcchimeravsLength) + "%**"
    except Exception as e:
         chimeraReads = "Problem reading outputfile"
intASV = 1
totalAsvs=""
intAsvs=1
try:
    asv_file=snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+"/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt"
    tasvs = subprocess.run( ["cat " +  asv_file + " | wc -l"], stdout=subprocess.PIPE, shell=True)
    intAsvs = int(tasvs.stdout.decode('utf-8').strip())
    #print("Total OTUS" + str(intOtus))
    totalAsvs = "**" + str(intAsvs) + "**"
except Exception as e:
    totalAsvs = "**Problem reading outputfile**"

prcAssigned = 0.0
prcNotAssignedOtus="TBD"
assignedOtus=0
notAssignedOtus=0
try:
    aOtus = subprocess.run( ["cat " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt | cut -f2 | grep -w NA | wc -l"], stdout=subprocess.PIPE, shell=True)
    notAssignedOtus = int(aOtus.stdout.decode('utf-8').strip())
    #print("Not assigned OTUS" + str(notAssignedOtus))
    assignedOtus = (intAsvs - notAssignedOtus)
    prcAssigned = float(assignedOtus/intAsvs)*100

    prcAssignedAsvs = "**" + "{:.2f}".format(prcAssigned) + "%**"
except Exception as e:
    prcAssignedAsvs = "**Problem reading outputfile**"


intSingletons = 1;
totalSingletons =""
try:
    totS = subprocess.run( ["grep -v \"^#\" " +  snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/asvTable_noSingletons.txt" + " | wc -l"], stdout=subprocess.PIPE, shell=True)
    intSingletons = int(totS.stdout.decode('utf-8').strip())
    #print("Total OTUS" + str(intOtus))
    totalSingletons = "**" + str(intSingletons) + "**"
except Exception as e:
    totalSingletons = "**Problem reading outputfile**"


notAssignedSingleOtus = 0
assignedSingleOtus = 0
totalAssignedSingletons =""
try:
    sOtus = subprocess.run( ["cat " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta  |  grep '^>' | sed 's/>//' | grep -F -w -f - " +snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt | cut -f2 | grep -w NA | wc -l" ], stdout=subprocess.PIPE, shell=True)
    notAssignedSingleOtus = int(sOtus.stdout.decode('utf-8').strip())
#print("Not assigned OTUS" + str(notAssignedOtus))
    assignedSingleOtus = (intSingletons - notAssignedSingleOtus)
    totalAssignedSingletons = "**" + str(assignedSingleOtus) + "%**"
except Exception as e:
    totalAssignedSingletons = "**Problem reading outputfile**"
  
prcSingle = 0.0
prcSingleStr=""  
try:
    prcSingle=float(assignedSingleOtus/intSingletons)*100
    prcSingleStr = "**" + "{:.2f}".format(prcSingle) + "%**" 
except Exception as e:
    prcSingleStr="**Error parsing output**"


#include user description on the report
desc = snakemake.config["description"]
txtDescription = ""
if len(desc) > 0:
    txtDescription = "\n**User description:** "+desc+"\n"


################################################################################
#                       Sample distribution chart                              #
################################################################################
countTxt="Following the read counts: \n\n"
fileData = []
headers = []
data =[]
headers.append("File description")
headers.append("Location")
headers.append("#")
headers.append("(%)")
fileData.append(headers)
#combined
data.append("Demultiplexed reads")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/<LIBRARY>_data/demultiplexed/\*.fastq.gz")
data.append(str(intTotalReads))
data.append("100%")
fileData.append(data)
data=[]
#filtered
data.append("QA filtered & trimmed reads")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/<LIBRARY>_data/demultiplexed/filtered/\*.fastq.gz")
data.append(str(intFilteredReads))
data.append("{:.2f}".format(float(prcFiltered))+"%")
fileData.append(data)
data=[]

#fw denoised
data.append("Denoised FW reads")
data.append("*No intermediate file generated*")
data.append(str(intDenoisedFWReads))
data.append("{:.2f}".format(prcDenoisedFW)+"%")
fileData.append(data)
data=[]

#rv denoised
data.append("Denoised RV reads")
data.append("*NO intermediate file generated*")
data.append(str(intDenoisedRVReads))
data.append("{:.2f}".format(prcDenoisedRV)+"%")
fileData.append(data)
data=[]

#Merged
data.append("Merged and full denoised reads")
data.append("*No intermediate file generated*")
data.append(str(intmergedReads))
data.append("{:.2f}".format(prcmerged)+"%")
fileData.append(data)
data=[]

#LengthFiltered
data.append("Length filtered")
data.append("*No intermediate file generated*")
data.append(str(intlengthFReads))
data.append("{:.2f}".format(prclengthF)+"%")
fileData.append(data)
data=[]

if snakemake.config["dada2_asv"]["chimeras"] == "T":
    data.append("Chimera removed")
    data.append("*No intermediate file generated*")
    data.append(str(intchimeraReads))
    data.append("{:.2f}".format(prcchimera)+"%")
    fileData.append(data)
    data=[]

#asv
data.append("ASV table")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/asvTable.txt")
data.append(str(intAsvs))
data.append("{:.2f}".format(float((intAsvs/intTotalReads)*100))+"%")
fileData.append(data)
data=[]
#Taxonomy
data.append("Taxonomy assignation")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt")
data.append(str(assignedOtus))
data.append("{:.2f}".format(float((assignedOtus/intAsvs)*100))+"%")
fileData.append(data)
data=[]
#otus no singletons
data.append("ASV table (no singletons: a > " + str(snakemake.config["filterOtu"]["n"])+")")
data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/asvTable_noSingletons.txt")
data.append(str(intSingletons))
data.append("{:.2f}".format(float((intSingletons/intAsvs)*100))+"%")
fileData.append(data)
data=[]
#Assigned singletons
data.append("Assigned no singletons")
data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/asvTable_noSingletons.txt")
data.append(str(assignedSingleOtus))
try:
    data.append("{:.2f}".format(prcSingle)+"%")
except Exception as e:
    data.append("Err")
    print("Error - Assigned no singletons - dividing: "+ str(assignedSingleOtus)+"/"+ str(intSingletons))
fileData.append(data)
countTxt += make_table(fileData)
################################################################################
#                         Generate sequence amounts chart                      #
################################################################################
numbers=[intTotalReads];
labels=["Initial\nreads"];
prcs=[]

prcs.append("100%")
#if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
#    numbers.append(intDerep)
#    labels.append("Derep.")
#    prcs.append("{:.2f}".format(float((intDerep/intTotalReads)*100))+"%")

numbers.append(intFilteredReads)
labels.append("Filtered\nreads")
prcs.append("{:.2f}".format(prcFiltered)+"%")

#numbers.append(intDenoisedFWReads)
#labels.append("Denoised\nFW reads")
#prcs.append("{:.2f}".format(prcDenoisedFW)+"%")

#numbers.append(intDenoisedRVReads)
#labels.append("Denoised\nRV reads")
#prcs.append("{:.2f}".format(prcDenoisedRV)+"%")


numbers.append(intmergedReads)
labels.append("Merged\nreads")
prcs.append("{:.2f}".format(prcmerged)+"%")

numbers.append(intlengthFReads)
labels.append("Length\nfiltered")
prcs.append("{:.2f}".format(prclengthF)+"%")

if snakemake.config["dada2_asv"]["chimeras"] == "T":
    numbers.append(intchimeraReads)
    labels.append("Chimera\nremoved")
    prcs.append("{:.2f}".format(prcchimera)+"%")

numbers.append(intAsvs)
labels.append("ASVs")
prcs.append("{:.2f}".format(float((intAsvs/intTotalReads)*100))+"%")

numbers.append(assignedOtus)
labels.append("Assigned\nASVs")
prcs.append("{:.2f}".format(float((assignedOtus/intAsvs)*100))+"%")

numbers.append(intSingletons)
labels.append("No\nSingletons")
prcs.append("{:.2f}".format(float((intSingletons/intAsvs)*100))+"%")

numbers.append(assignedSingleOtus)
labels.append("Assigned no\nsingletons")
prcs.append("{:.2f}".format(prcSingle)+"%")

createChartPrc(numbers, tuple(labels),prcs,snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files/sequence_numbers_asv.png")

###############################################################################
#                       Varaible sections                                     #
################################################################################
variable_refs=""
assignTaxoStr = ""
if snakemake.config["ANALYSIS_TYPE"] == "ASV":
    assignTaxoStr =":red:`Tool:` RDP_\n\n"
    assignTaxoStr += ":green:`Function:` assignTaxonomy() *implementation of RDP Classifier within dada2*\n\n"
    assignTaxoStr += ":green:`Reference database:` " + str(snakemake.config["dada2_taxonomy"]["db"])+ "\n\n"
    if snakemake.config["dada2_taxonomy"]["add_sps"]["add"].casefold() == "T":
        assignTaxoStr += ":green:`Taxonomy species file:` " + str(snakemake.config["dada2_taxonomy"]["add_sps"]["db_sps"])+ "\n\n"
    variable_refs+=".. [RDP]  Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.\n\n"
elif snakemake.config["assignTaxonomy"]["tool"] == "blast":
    assignTaxoStr =":red:`Tool:` RDP_\n\n"
    assignTaxoStr += ":red:`Version:` " + blastnVersion + "\n\n"
    variable_refs+= ".. [blast] Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215(3):403-410\n\n"
    ref = ""
    if len(str(snakemake.config["assignTaxonomy"]["blast"]["blast_db"])) > 1:
        assignTaxoStr +=  ":green:`Reference database:` "+ str(snakemake.config["assignTaxonomy"]["blast"]["blast_db"])+"\n\n"
        ref= "-db " + str(snakemake.config["assignTaxonomy"]["blast"]["blast_db"])
    else:
        assignTaxoStr +=  ":green:`Reference fasta file:` "+ str(snakemake.config["assignTaxonomy"]["blast"]["fasta_db"])+"\n\n"
        ref= "-subject "+ str(snakemake.config["assignTaxonomy"]["blast"]["fasta_db"])
    assignTaxoStr +=  ":green:`Taxonomy mapping file:` "+ str(snakemake.config["assignTaxonomy"]["blast"]["mapFile"])+"\n\n"
    assignTaxoStr += "**Command:**\n\n"
    assignTaxoStr += ":commd:`"+ str(snakemake.config["assignTaxonomy"]["blast"]["command"] )+" " +ref + "-evalue " + str(snakemake.config["assignTaxonomy"]["blast"]["evalue"]) + "-outfmt '6 qseqid sseqid pident qcovs evalue bitscore' -num_threads " + str(snakemake.config["assignTaxonomy"]["blast"]["jobs"]) + " -max_target_seqs "
    assignTaxoStr += str(snakemake.config["assignTaxonomy"]["blast"]["max_target_seqs"]) +" -perc_identity "+ str(snakemake.config["assignTaxonomy"]["blast"]["identity"]) + " -out representative_seq_set_tax_blastn.out`\n\n"
    if snakemake.config["assignTaxonomy"]["blast"]["max_target_seqs"] != 1:
        assignTaxoStr += "After blast assignation, **results were mapped to their LCA using stampa_merge.py** script\n\n"

elif snakemake.config["assignTaxonomy"]["tool"] == "qiime":
    assignTaxoStr =":red:`Tool:` [QIIME]_\n\n"
    assignTaxoStr += ":red:`Version:` "+assignTaxaVersion
    assignTaxoStr += ":green:`Method:` **" + str(snakemake.config["assignTaxonomy"]["qiime"]["method"])+ "**\n\n"
    assignTaxoStr += "Reference database: " + str(snakemake.config["assignTaxonomy"]["qiime"]["dbFile"])+ "\n\n"
    assignTaxoStr += "Taxonomy mapping file: " + str(snakemake.config["assignTaxonomy"]["qiime"]["mapFile"])+ "\n\n"
    assignTaxoStr += "**Command:**\n\n"
    assignTaxoStr += ":commd:`parallel_assign_taxonomy_" + str(snakemake.config["assignTaxonomy"]["qiime"]["method"])+ ".py -i " + str(snakemake.wildcards.PROJECT)+ "/runs/" + str(snakemake.wildcards.run)+ "/otu/representative_seq_set.fasta --id_to_taxonomy_fp " + str(snakemake.config["assignTaxonomy"]["qiime"]["mapFile"])+ " --reference_seqs_fp "
    assignTaxoStr += str(snakemake.config["assignTaxonomy"]["qiime"]["dbFile"])+ " --jobs_to_start " + str(snakemake.config["assignTaxonomy"]["qiime"]["jobs"])+ " " + str(snakemake.config["assignTaxonomy"]["qiime"]["extra_params"])+ " "
    assignTaxoStr += "--output_dir " + str(snakemake.wildcards.PROJECT)+ "/runs/" + str(snakemake.wildcards.run)+ "/otu/taxonomy_" + str(snakemake.config["assignTaxonomy"]["tool"])+ "/`\n\n"
elif snakemake.config["assignTaxonomy"]["tool"] == "vsearch":
    assignTaxoStr =":red:`Tool:` [vsearch]_\n\n"
    assignTaxoStr += ":red:`Version:` " + vsearchVersion_tax + "\n\n"
    assignTaxoStr +=  ":green:`Reference fasta file:` "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["db_file"])+"\n\n"
    assignTaxoStr +=  ":green:`Taxonomy mapping file:` "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["mapFile"])+"\n\n"
    assignTaxoStr += "**Command:**\n\n"
    assignTaxoStr += ":commd:`"+ str(snakemake.config["assignTaxonomy"]["vsearch"]["command"] )+ "--usearch_global "+ str(snakemake.wildcards.PROJECT)+ "/runs/" + str(snakemake.wildcards.run)+ "/otu/representative_seq_set.fasta --db "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["db_file"])
    assignTaxoStr += " --dbmask none --qmask none --rowlen 0 --id "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["identity"])+" --iddef " + str(snakemake.config["assignTaxonomy"]["vsearch"]["identity_definition"])+" --userfields query+id" + str(snakemake.config["assignTaxonomy"]["vsearch"]["identity_definition"])+"+target "
    assignTaxoStr += " --maxaccepts "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["max_target_seqs"]) + " --threads " + str(snakemake.config["assignTaxonomy"]["vsearch"]["jobs"]) + " "+ str(snakemake.config["assignTaxonomy"]["vsearch"]["extra_params"]) + " --output_no_hits --userout  representative_seq_set_tax_vsearch.out`\n\n"
    if (snakemake.config["assignTaxonomy"]["vsearch"]["max_target_seqs"]) != 1:
        assignTaxoStr += "After vsearch assignation, **results were mapped to their LCA using stampa_merge.py** script\n\n"

#Alignment report
alignmentReport = ""
if snakemake.config["alignRep"]["align"] == "T":
    alignmentReport = "\nAlign representative sequences\n-------------------------------\n\n"
    alignmentReport+="Align the sequences in a FASTA file to each other or to a template sequence alignment.\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - align_seqs.py\n\n"
    alignmentReport+=":red:`Version:` "+alignFastaVersion +"\n\n"
    alignmentReport+=":green:`Method:` ["+ snakemake.config["alignRep"]["m"] + "]_\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`align_seqs.py -m "+snakemake.config["alignRep"]["m"] +" -i "+ snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2/representative_seq_set_noSingletons.fasta "+ snakemake.config["alignRep"]["extra_params"] + " -o "
    alignmentReport+=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta`\n\n"
    alignmentReport+="**Output files:**\n\n"
    alignmentReport+=":green:`- Aligned fasta file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta\n\n"
    alignmentReport+=":green:`- Log file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_log.txt\n\n"
    alignmentReport+=alignBenchmark+"\n\n"

    alignmentReport+="Filter alignment\n-----------------\n\n"
    alignmentReport+="Removes positions which are gaps in every sequence.\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - filter_alignment.py\n\n"
    alignmentReport+=":red:`Version:` "+filterAlignmentVersion +"\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`filter_alignment.py -i  "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta " +snakemake.config["filterAlignment"]["extra_params"]
    alignmentReport+=" -o  "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/filtered/`\n\n"
    alignmentReport+="**Output file:**\n\n"
    alignmentReport+=":green:`- Aligned fasta file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned_pfiltered.fasta\n\n"
    alignmentReport+=alignFilteredBenchmark+"\n\n"

    alignmentReport+="Make tree\n-----------\n\n"
    alignmentReport+="Create phylogenetic tree (newick format).\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - make_phylogeny.py\n\n"
    alignmentReport+=":red:`Version:` "+makePhyloVersion +"\n\n"
    alignmentReport+=":green:`Method:` ["+ snakemake.config["makeTree"]["method"] + "]_\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`make_phylogeny.py -i "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.fasta -o representative_seq_set_noSingletons_aligned_pfiltered.tre "+ snakemake.config["makeTree"]["extra_params"]+ " -t " + snakemake.config["makeTree"]["method"]+"`\n\n"
    alignmentReport+="**Output file:**\n\n"
    alignmentReport+=":green:`- Taxonomy tree:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/aligned/representative_seq_set_noSingletons_aligned.tre\n\n"
    alignmentReport+=makePhyloBenchmark+"\n\n"
#KRONA REPORT
kronaReport = ""
if  snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kronaReport+="Krona report\n----------------\n\n"
    kronaReport+="Krona allows hierarchical data to be explored with zooming, multi-layered pie charts.\n\n"
    kronaReport+=":red:`Tool:` [Krona]_\n\n"
    if snakemake.config["krona"]["otu_table"].casefold() != "singletons":
        kronaReport+="These charts were created using the ASV table **without** singletons\n\n"
    else:
        kronaReport+="These charts were created using the ASV table **including** singletons\n\n"

    if snakemake.config["krona"]["samples"].strip() == "all":
        kronaReport+="The report was executed for all the samples.\n\n"
    else:
        kronaReport+="The report was executed for the following target samples: "+ snakemake.config["krona"]["samples"].strip() + "\n\n"
    if "-c" in snakemake.config["krona"]["extra_params"]:
        kronaReport+="The samples were combined on a single chart\n\n"
    else:
        kronaReport+="Each sample is represented on a separated chart (same html report).\n\n"
    kronaReport+="You can see the report at the following link:\n\n"
    kronaReport+=":green:`- Krona report:` kreport_\n\n"
    #kronaReport+=" .. _kreport: ../../runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/krona_report.html\n\n"
    kronaReport+=" .. _kreport: report_files/krona_report.dada2.html\n\n"

    kronaReport+="Or access the html file at:\n\n"
    kronaReport+=":green:`- Krona html file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/taxonomy_dada2/krona_report.html\n\n"
    kronaReport+=kronaBenchmark+"\n\n"

###############################################################################
#                         REFERENCES                                     #
################################################################################
#dada2
variable_refs+= ".. [dada2] Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13, 581-583. doi: 10.1038/nmeth.3869.\n\n"

#ALIGNMENT
if snakemake.config["alignRep"]["align"] == "T":
    if snakemake.config["alignRep"]["m"] == "pynast":
        variable_refs+= ".. [pynast] Caporaso JG, Bittinger K, Bushman FD, DeSantis TZ, Andersen GL, Knight R. 2010. PyNAST: a flexible tool for aligning sequences to a template alignment. Bioinformatics 26:266-267.\n\n"
    elif snakemake.config["alignRep"]["m"] == "infernal":
        variable_refs+= ".. [infernal] Nawrocki EP, Kolbe DL, Eddy SR. 2009. Infernal 1.0: Inference of RNA alignments. Bioinformatics 25:1335-1337.\n\n"

    if snakemake.config["makeTree"]["method"] == "fasttree":
        variable_refs+= ".. [fasttree] Price MN, Dehal PS, Arkin AP. 2010. FastTree 2-Approximately Maximum-Likelihood Trees for Large Alignments. Plos One 5(3).\n\n"
    elif snakemake.config["makeTree"]["method"] == "raxml":
        variable_refs+= "..[raxml] Stamatakis A. 2006. RAxML-VI-HPC: Maximum Likelihood-based Phylogenetic Analyses with Thousands of Taxa and Mixed Models. Bioinformatics 22(21):2688-2690.\n\n"
    elif snakemake.config["makeTree"]["method"] == "clearcut":
        variable_refs+= "..[clearcut] Evans J, Sheneman L, Foster JA. 2006. Relaxed Neighbor-Joining: A Fast Distance-Based Phylogenetic Tree Construction Method. J Mol Evol 62:785-792.\n\n"
    elif snakemake.config["makeTree"]["method"] == "clustalw":
        variable_refs+= "..[clustalw] Larkin MA, Blackshields G, Brown NP, Chenna R, McGettigan PA, McWilliam H, Valentin F, Wallace IM, Wilm A, Lopez R, Thompson JD, Gibson TJ, Higgins DG. 2007. Clustal W and Clustal X version 2.0. Bioinformatics 23:2947-2948.\n\n"

########
# EXTRA
##############

errorPlots="" 
if snakemake.config["dada2_asv" ]["generateErrPlots"].casefold() == "t" or snakemake.config["dada2_asv" ]["generateErrPlots"].casefold() == "true":
    errorPlots+="**Error plots:** \n\n:green:`- FW reads error plot::`  " + snakemake.wildcards.PROJECT + "/runs/"+snakemake.wildcards.run+ "/asv/fw_err.pdf\n\n" 
    errorPlots+=":green:`- RV reads error plot::`  " + snakemake.wildcards.PROJECT + "/runs/"+snakemake.wildcards.run+ "/asv/rv_err.pdf\n\n"

#shorts and longs
shorts = str(snakemake.config["rm_reads"]["shorts"])
longs = str(snakemake.config["rm_reads"]["longs"])
with open(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/shorts_longs.log") as trimlog:
    i=0
    for line in trimlog:
        i=i+1
        #tokens = line.split("\t")
        if i== 1:
            shorts = line
        else:
            longs = line

chimeras="" 
if snakemake.config["dada2_asv" ]["chimeras"].casefold() == "t" or snakemake.config["dada2_asv" ]["chimeras"].casefold() == "true":
    chimeras="Remove chimeras\n~~~~~~~~~~~~~~~~\n\n"
    chimeras+="Sequence variants identified as bimeric are removed, and a bimera-free collection of unique sequences is generated.\n\n"
    chimeras+=":green:`Function:` removeBimeraDenovo()\n\n"
    chimeras+=":green:`Method:` consensus\n\n" 

report("""
{title}
    .. role:: commd
    .. role:: red
    .. role:: green

**CASCABEL** is designed to run amplicon sequence analysis across single or multiple read libraries. This report consists of the ASV table creation and taxonomic assignment for all the combined accepted reads of given samples or libraries, if multiple.

{txtDescription}

Filter and Trim
---------------
Once that all the individual libraries were demultiplexed, the fastq files from all the samples for all the libraries were processed together. 

The filter and trimming steps were both performed with the **filterAndTrim()** function from the R package dada2, according to user parameters.

:red:`Tool:` dada2_ 

:red:`Version:` {dada2Version}

:green:`Function:` filterAndTrim()

:green:`Max Expected Errors (maxEE) FW:` {snakemake.config[dada2_filter][maxEE_FW]}

:green:`Max Expected Errors (maxEE) RV:` {snakemake.config[dada2_filter][maxEE_RV]}

:green:`Forward read truncation:` {snakemake.config[dada2_filter][truncFW]}

:green:`Reverse read truncation:` {snakemake.config[dada2_filter][truncRV]}

**Command:**


:commd:`Scripts/asvFilter.R $PWD {snakemake.config[dada2_filter][generateQAplots]} {snakemake.config[dada2_filter][truncFW]} {snakemake.config[dada2_filter][truncRV]} {snakemake.config[dada2_filter][maxEE_FW]} {snakemake.config[dada2_filter][maxEE_RV]} {snakemake.config[dada2_filter][cpus]} {snakemake.config[dada2_filter][extra_params]} {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/filter_summary.out`


**Output file:**

:green:`- Filtered fastq files:`   {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/<Library>/demultiplexed/filtered/

:green:`- Summary:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/filter_summary.out


:red:`Note:` To speed up downstream computation, consider tightening maxEE. If too few reads are passing the filter, consider relaxing maxEE, perhaps especially on the reverse reads.

Make sure that your forward and reverse reads overlap after length truncation.

{asvFilterBenchmark}

 
Amplicon Sequence Variants
----------------------------
In order to identify ASVs, dada2 workflow require to execute several steps. Following a summary of these steps and its main parameters. 

:red:`Tool:` dada2_ 

:red:`Version:` {dada2Version}

Learn errors
~~~~~~~~~~~~~~~~
The first step after filtering the reads is to learn the errors from the fastq files.

:green:`Function:` learnErrors(filteredFQ)

{errorPlots}

ASV inference
~~~~~~~~~~~~~~~
The amplicon sequence variant identification consists of a high resolution sample inference from the amplicon data using the learned errors. 
 
:green:`Function:` dada(filteredFQ, errors, pool='{snakemake.config[dada2_asv][pool]}')
 
Merge pairs
~~~~~~~~~~~~~~~
In this step, forward and reverse reads are paired in order to create full denoised sequences.

:green:`Function:` mergePairs(dadaF, dadaR)

:green:`Min overlap:` {snakemake.config[dada2_merge][minOverlap]}

:green:`Max mismatch:` {snakemake.config[dada2_merge][maxMismatch]}

Length filtering   
~~~~~~~~~~~~~~~~~~
Sequences that are much longer or shorter than expected may be the result of non-specific priming.

:green:`- Shortest length:` {shorts}

:green:`- Longest length:` {longs}

{chimeras}

**Output files:**

:green:`- Representative ASV sequences:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/representative_seq_set.fasta

The total number of different ASVs is: {totalAsvs}


Assign taxonomy
----------------
Given a set of sequences, assign the taxonomy of each sequence.

{assignTaxoStr}

The percentage of successfully assigned ASVs is: {prcAssignedAsvs}

**Output file:**

:green:`- ASV taxonomy assignation:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt


The previous steps were performed within a Cascabel R script according to the following command:

**Command**

:commd:`Scripts/asvDada2.R $PWD  {snakemake.config[dada2_asv][pool]}   {snakemake.config[dada2_asv][cpus]}    {snakemake.config[dada2_asv][generateErrPlots]}   {snakemake.config[dada2_asv][extra_params]}  {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/    {snakemake.config[rm_reads][shorts]}    {snakemake.config[rm_reads][longs]}   {snakemake.config[rm_reads][offset]}    {snakemake.config[dada2_asv][chimeras]}    {snakemake.config[dada2_taxonomy][db]}   {snakemake.config[dada2_taxonomy][add_sps][db_sps]}    {snakemake.config[dada2_taxonomy][add_sps][add]}   {snakemake.config[dada2_taxonomy][extra_params]}  {snakemake.config[dada2_merge][minOverlap]}  {snakemake.config[dada2_merge][maxMismatch]}  {snakemake.config[dada2_taxonomy][add_sps][extra_params]}`  


{dada2Benchmark}

Make ASV table
---------------
Tabulates the number of times an ASV is found in each sample, and adds the taxonomic predictions for each ASV in the last column.

**Command:**

:commd:`cat {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt | awk 'NR==FNR{{if(NR>1){{tax=$2;for(i=3;i<=NF;i++){{tax=tax";"$i}};h[$1]=tax;}}next;}} {{if(FNR==1){{print $0"\\ttaxonomy"}}else{{print $0"\\t"h[$1]}}' - {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/asv_table.txt  >  {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.txt`

**Output file:**

:green:`- ASV table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.txt

{otuTableBenchmark}

Convert ASV table
------------------
Convert from txt to the BIOM table format.

:red:`Tool:` [BIOM]_

:red:`Version:` {convertBiomVersion}

**Command:**

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.txt -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.biom {snakemake.config[biom][tableType]} --table type "OTU table"  --to-hdf5 --process-obs-metdata taxonomy`

**Output file:**

:green:`- Biom format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.biom

{convertOtuBenchmark}

Summarize Taxa
---------------
Summarize information of the representation of taxonomic groups within each sample.

:red:`Tool:` [QIIME]_ - summarize_taxa.py

:red:`Version:` {summTaxaVersion}

**Command:**

:commd:`summarize_taxa.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/otuTable.biom {snakemake.config[summTaxa][extra_params]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/summary/`

**Output file:**

:green:`- Taxonomy summarized counts at different taxonomy levels:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/summary/otuTable_L**N**.txt

Where **N** is the taxonomy level. Default configuration produces levels from 2 to 6.

{summTaxaBenchmark}

Filter ASV table
-----------------
Filter ASVs from an ASV table based on their observed counts or identifier.

:red:`Tool:` [QIIME]_ - filter_otus_from_otu_table.py

:red:`Version:` {filterOTUNoSVersion}

:green:`Minimum observation counts:` {snakemake.config[filterOtu][n]}

**Command:**

:commd:`filter_otus_from_otu_table.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable_noSingletons.biom {snakemake.config[filterOtu][extra_params]} -n {snakemake.config[filterOtu][n]}`

**Output file:**

:green:`- Biom table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/otuTable_noSingletons.biom

{asvNoSingletonsBenchmark}

Convert Filtered ASV table
---------------------------
Convert the filtered OTU table from the BIOM table format to a human readable format

:red:`Tool:` [BIOM]_

:red:`Version:` {convertBiomVersion}

**Command:**

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_dada2/asvTable_noSingletons.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable_noSingletons.txt {snakemake.config[biom][tableType]} {snakemake.config[biom][headerKey]} {snakemake.config[biom][extra_params]} {snakemake.config[biom][outFormat]}`

**Output file:**

:green:`- TSV format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/asv/taxonomy_dada2/asvTable_noSingletons.txt

{filterASVTableBenchmark}

Filter representative sequences
---------------------------------
Remove sequences according to the filtered OTU biom table.

:red:`Tool:` [QIIME]_ - filter_fasta.py

:red:`Version:` {filterFastaVersion}

**Command:**

:commd:`filter_fasta.py -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/asv/representative_seq_set.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta {snakemake.config[filterFasta][extra_params]} -b {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/asv/taxonomy_dada2/otuTable_noSingletons.biom`

**Output file:**

:green:`- Filtered fasta file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/asv/taxonomy_dada2/representative_seq_set_noSingletons.fasta


{alignmentReport}

{kronaReport}

Final counts
-------------

{countTxt}

.. image:: report_files/sequence_numbers_asv.png

:red:`Note:`

:green:`- Assigned ASVs percentage` is the amount of successfully assigned ASVs.

:green:`- No singletons percentage` is the percentage of no singletons ASVs in reference to the complete ASV table.

:green:`- Assigned No singletons` is the amount of successfully no singletons assigned ASVs.

References
------------

.. [QIIME] QIIME. Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Gonzalez Pena A, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. 2010. QIIME allows analysis of high-throughput community sequencing data. Nature Methods 7(5): 335-336.

.. [Cutadapt] Cutadapt v1.15 .Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

.. [vsearch] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

.. [Krona] Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

.. [BIOM] The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso.GigaScience 2012, 1:7. doi:10.1186/2047-217X-1-7

{variable_refs}


""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
