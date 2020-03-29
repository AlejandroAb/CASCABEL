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


title = getSampleList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/samples.log")
catCommand =  getCombinedSamplesList(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/cat_samples.log")

################################################################################
#                         Benchmark Section                                    #
# This section is to generate a pre-formatted text with the benchmark info for #
# All the executed rules.                                                      #
################################################################################
combineBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/combine_seqs_fw_rev.benchmark")
otuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu.benchmark")
pikRepBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/pick_reps.benchmark")
assignTaxaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/assign_taxa.benchmark")
otuTableBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable.biom.benchmark")
convertOtuBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable.txt.benchmark")
summTaxaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/summary/summarize_taxa.benchmark")
otuNoSingletonsBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable_nosingletons.bio.benchmark")
filterBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/representative_seq_set_noSingletons.benchmark")
deRepBenchmark=""
if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
    deRepBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/derep.benchmark")
if snakemake.config["alignRep"]["align"] == "T":
    #align_seqs.py -m {config[alignRep][m]} -i {input} -o {params.outdir} {config[alignRep][extra_params]}
    alignBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/align_rep_seqs.benchmark")
    #"filter_alignment.py -i {input} -o {params.outdir} {config[filterAlignment][extra_params]}"
    alignFilteredBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/filtered/align_rep_seqs.benchmark")
    #"make_phylogeny.py -i {input} -o {output} {config[makeTree][extra_params]}"
    makePhyloBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/filtered/representative_seq_set_noSingletons_aligned_pfiltered.benchmark")
kronaBenchmark=""
if snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kronaBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/krona_report.benchmark")

#dada2FilterBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/filter.benchmark")
#dada2Benchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2.benchmark")
#dada2BiomBenchmark = readBenchmark(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/asv/dada2.biom.benchmark")

################################################################################
#                         TOOLS VERSION SECTION                          #
################################################################################
clusterOtuV = subprocess.run([snakemake.config["qiime"]["path"]+'pick_otus.py', '--version'], stdout=subprocess.PIPE)
clusterOtuVersion = "**" + clusterOtuV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

pickRepV = subprocess.run([snakemake.config["qiime"]["path"]+'pick_rep_set.py', '--version'], stdout=subprocess.PIPE)
pickRepVersion = "**" + pickRepV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

assignTaxaV = subprocess.run([snakemake.config["qiime"]["path"]+'parallel_assign_taxonomy_'+snakemake.config["assignTaxonomy"]["qiime"]["method"]+'.py', '--version'], stdout=subprocess.PIPE)
assignTaxaVersion = "**" + assignTaxaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

makeOTUV = subprocess.run([snakemake.config["qiime"]["path"]+'make_otu_table.py', '--version'], stdout=subprocess.PIPE)
makeOTUVersion = "**" + makeOTUV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

convertBiomV = subprocess.run([snakemake.config["biom"]["command"], '--version'], stdout=subprocess.PIPE)
convertBiomVersion = "**" + convertBiomV.stdout.decode('utf-8').strip() + "**"

summTaxaSV = subprocess.run([snakemake.config["qiime"]["path"]+'summarize_taxa.py', '--version'], stdout=subprocess.PIPE)
summTaxaVersion = "**" + summTaxaSV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

filterOTUNoSV = subprocess.run([snakemake.config["qiime"]["path"]+'filter_otus_from_otu_table.py', '--version'], stdout=subprocess.PIPE)
filterOTUNoSVersion = "**" + filterOTUNoSV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

filterFastaV = subprocess.run([snakemake.config["qiime"]["path"]+'filter_fasta.py', '--version'], stdout=subprocess.PIPE)
filterFastaVersion = "**" + filterFastaV.stdout.decode('utf-8').replace('Version:','').strip() + "**"

blastnV = subprocess.run([snakemake.config["assignTaxonomy"]["blast"]["command"], '-version'], stdout=subprocess.PIPE)
blastnVersion = "**" + blastnV.stdout.decode('utf-8').split('\n', 1)[0].replace('blastn:','').strip() + "**"

vsearchV2 = subprocess.run([snakemake.config["assignTaxonomy"]["vsearch"]["command"], '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
vsearchVersion_tax = "**" + vsearchV2.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"

if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
    vsearchV = subprocess.run([snakemake.config["derep"]["vsearch_cmd"], '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    vsearchVersion = "**" + vsearchV.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"

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
    treads = subprocess.run( ["grep '^>' " + snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/seqs_fw_rev_combined.fasta | wc -l"], stdout=subprocess.PIPE, shell=True)
    intTotalReads = int(treads.stdout.decode('utf-8').strip())
    totalReads = "**" + str(intTotalReads) + "**"
except Exception as e:
    totalReads = "Problem reading outputfile"

derep_reads = "TBD"
intDerep=1
if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
    try:
        totd = subprocess.run( ["grep \"^>\" " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/derep/seqs_fw_rev_combined_derep.fasta" + " | wc -l"], stdout=subprocess.PIPE, shell=True)
        intDerep = int(totd.stdout.decode('utf-8').strip())
        derep_reads = "**" + str(intDerep) + "**"
    except Exception as e:
        derep_reads = "**Problem reading outputfile**"

intOtus = 1
try:
    otu_file=""
    if (snakemake.config["derep"]["dereplicate"] == "T" and snakemake.config["pickOTU"]["m"] != "swarm" and snakemake.config["pickOTU"]["m"] != "usearch"):
        otu_file = snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/seqs_fw_rev_combined_remapped_otus.txt"
    else:
        otu_file = snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/seqs_fw_rev_combined_otus.txt"
    totus = subprocess.run( ["cat " +  otu_file + " | wc -l"], stdout=subprocess.PIPE, shell=True)
    intOtus = int(totus.stdout.decode('utf-8').strip())
    #print("Total OTUS" + str(intOtus))
    totalOtus = "**" + str(intOtus) + "**"
except Exception as e:
    totalOtus = "**Problem reading outputfile**"

prcAssigned = 0.0
prcNotAssignedOtus="TBD"
try:
    nohit = "'No blast hit|Unassigned'"
    #if snakemake.config["assignTaxonomy"]["tool"] != "blast":
    #    nohit = "'Unassigned'"
    aOtus = subprocess.run( ["grep -E "+ nohit + " " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/representative_seq_set_tax_assignments.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    notAssignedOtus = int(aOtus.stdout.decode('utf-8').strip())
    #print("Not assigned OTUS" + str(notAssignedOtus))
    assignedOtus = (intOtus - notAssignedOtus)
    prcAssigned = (assignedOtus/intOtus)*100

    prcAssignedOtus = "**" + "{:.2f}".format(prcAssigned) + "%**"
except Exception as e:
    prcAssignedOtus = "**Problem reading outputfile**"


intSingletons = 1;
try:
    totS = subprocess.run( ["grep -v \"^#\" " +  snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable_noSingletons.txt" + " | wc -l"], stdout=subprocess.PIPE, shell=True)
    intSingletons = int(totS.stdout.decode('utf-8').strip())
    #print("Total OTUS" + str(intOtus))
    totalSingletons = "**" + str(intSingletons) + "**"
except Exception as e:
    totalSingletons = "**Problem reading outputfile**"

nohit = "'No blast hit|Unassigned|None'"
#if snakemake.config["assignTaxonomy"]["tool"] != "blast":
#    nohit = "'Unassigned'"
notAssignedSingleOtus = 0
assignedSingleOtus = 0
try:
    sOtus = subprocess.run( ["grep -E "+ nohit + " " +  snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable_noSingletons.txt | wc -l"], stdout=subprocess.PIPE, shell=True)
    notAssignedSingleOtus = int(sOtus.stdout.decode('utf-8').strip())
#print("Not assigned OTUS" + str(notAssignedOtus))
    assignedSingleOtus = (intSingletons - notAssignedSingleOtus)
except Exception as e:
    totalAssignedSingletons = "**Problem reading outputfile**"


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
data.append("Combined clean reads")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/seqs_fw_rev_combined.fasta")
data.append(str(intTotalReads))
data.append("100%")
fileData.append(data)
data=[]
#derep
if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
	data.append("Dereplicated reads")
	data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/derep/seqs_fw_rev_combined_derep.fasta")
	data.append(str(intDerep))
	data.append("{:.2f}".format(float((intDerep/intTotalReads)*100))+"%")
	fileData.append(data)
	data=[]

#otus
data.append("OTU table")
data.append(otu_file)
data.append(str(intOtus))
data.append("{:.2f}".format(float((intOtus/intTotalReads)*100))+"%")
fileData.append(data)
data=[]
#Taxonomy
data.append("Taxonomy assignation")
data.append(snakemake.wildcards.PROJECT+ "/runs/" + snakemake.wildcards.run+ "/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/representative_seq_set_tax_assignments.txt")
data.append(str(assignedOtus))
data.append("{:.2f}".format(float((assignedOtus/intOtus)*100))+"%")
fileData.append(data)
data=[]
#otus no singletons
data.append("OTU table (no singletons: a > " + str(snakemake.config["filterOtu"]["n"])+")")
data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable_noSingletons.txt")
data.append(str(intSingletons))
data.append("{:.2f}".format(float((intSingletons/intOtus)*100))+"%")
fileData.append(data)
data=[]
#Assigned singletons
data.append("Assigned no singletons")
data.append(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/otuTable_noSingletons.txt")
data.append(str(assignedSingleOtus))
try:
    data.append("{:.2f}".format(float((assignedSingleOtus/intSingletons)*100))+"%")
except Exception as e:
    data.append("Err")
    print("Error - Assigned no singletons - dividing: "+ str(assignedSingleOtus)+"/"+ str(intSingletons))
fileData.append(data)
countTxt += make_table(fileData)
################################################################################
#                         Generate sequence amounts chart                      #
################################################################################
numbers=[intTotalReads];
labels=["Combined\nreads"];
prcs=[]

prcs.append("100%")
if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
    numbers.append(intDerep)
    labels.append("Derep.")
    prcs.append("{:.2f}".format(float((intDerep/intTotalReads)*100))+"%")

numbers.append(intOtus)
labels.append("OTUs")
prcs.append("{:.2f}".format(float((intOtus/intTotalReads)*100))+"%")

numbers.append(assignedOtus)
labels.append("Assigned\nOTUs")
prcs.append("{:.2f}".format(float((assignedOtus/intOtus)*100))+"%")

numbers.append(intSingletons)
labels.append("No\nSingletons")
prcs.append("{:.2f}".format(float((intSingletons/intOtus)*100))+"%")

numbers.append(assignedSingleOtus)
labels.append("Assigned NO\n singletons")
prcs.append("{:.2f}".format(float((assignedSingleOtus/intSingletons)*100))+"%")

createChartPrc(numbers, tuple(labels),prcs,snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/report_files/sequence_numbers_all.png")

###############################################################################
#                       Varaible sections                                     #
################################################################################
variable_refs=""
assignTaxoStr = ""
if snakemake.config["assignTaxonomy"]["tool"] == "blast":
    assignTaxoStr =":red:`Tool:` ["+str(snakemake.config["assignTaxonomy"]["tool"])+"]_\n\n"
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
        assignTaxoStr += "After taxonomy assignation with vsearch, top hits with the same sequence identity but different taxonomy were mapped to their last common ancestor (LCA) using the script **stampa_merge.py** from https://github.com/frederic-mahe/stampa.\n\n"

#Dereplication report
dereplicateReport=""
if  snakemake.config["derep"]["dereplicate"] == "T" and  snakemake.config["pickOTU"]["m"] != "swarm" and  snakemake.config["pickOTU"]["m"] != "usearch":
    dereplicateReport="Dereplicate reads\n"
    dereplicateReport+="---------------------\n\n"
    dereplicateReport+="Clusterize the reads with an identity threshold of 100%.\n\n"
    dereplicateReport+=":red:`Tool:` [vsearch]_\n\n"
    dereplicateReport+=":red:`Version:` " + vsearchVersion+"\n\n"
    dereplicateReport+="**Command:**\n\n"
    dereplicateReport+=":commd:`"+str(snakemake.config["derep"]["vsearch_cmd"]) +" --derep_fulllength  seqs_fw_rev_combined.fasta --output seqs_fw_rev_combined_derep.fasta --uc  seqs_fw_rev_combined_derep.uc --strand " + str(snakemake.config["derep"]["strand"]) + " --fasta_width 0 --minuniquesize "+ str(snakemake.config["derep"]["min_abundance"])+"`\n\n"
    dereplicateReport+="**Output files:**\n\n"
    dereplicateReport+=":green:`- Dereplicated fasta file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/seqs_fw_rev_combined_derep.fasta\n\n"
    dereplicateReport+=":green:`- Cluster file:` "+ snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/derep/seqs_fw_rev_combined_derep.uc\n\n"
    dereplicateReport+="Total number of dereplicated sequences is: "+str(derep_reads).strip()+"\n\n"+deRepBenchmark+"\n\n"
#Alignment report
alignmentReport = ""
if snakemake.config["alignRep"]["align"] == "T":
    alignmentReport = "\nAlign representative sequences\n-------------------------------\n\n"
    alignmentReport+="Align the sequences in a FASTA file to each other or to a template sequence alignment.\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - align_seqs.py\n\n"
    alignmentReport+=":red:`Version:` "+alignFastaVersion +"\n\n"
    alignmentReport+=":green:`Method:` ["+ snakemake.config["alignRep"]["m"] + "]_\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`align_seqs.py -m "+snakemake.config["alignRep"]["m"] +" -i "+ snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/"+snakemake.config["assignTaxonomy"]["tool"]+"/representative_seq_set_noSingletons.fasta "+ snakemake.config["alignRep"]["extra_params"] + " -o "
    alignmentReport+=snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned.fasta`\n\n"
    alignmentReport+="**Output files:**\n\n"
    alignmentReport+=":green:`- Aligned fasta file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned.fasta\n\n"
    alignmentReport+=":green:`- Log file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_log.txt\n\n"
    alignmentReport+=alignBenchmark+"\n\n"

    alignmentReport+="Filter alignment\n-----------------\n\n"
    alignmentReport+="Removes positions which are gaps in every sequence.\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - filter_alignment.py\n\n"
    alignmentReport+=":red:`Version:` "+filterAlignmentVersion +"\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`filter_alignment.py -i  "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned.fasta " +snakemake.config["filterAlignment"]["extra_params"]
    alignmentReport+=" -o  "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/filtered/`\n\n"
    alignmentReport+="**Output file:**\n\n"
    alignmentReport+=":green:`- Aligned fasta file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned_pfiltered.fasta\n\n"
    alignmentReport+=alignFilteredBenchmark+"\n\n"

    alignmentReport+="Make tree\n-----------\n\n"
    alignmentReport+="Create phylogenetic tree (newick format).\n\n"
    alignmentReport+=":red:`Tool:` [QIIME]_ - make_phylogeny.py\n\n"
    alignmentReport+=":red:`Version:` "+makePhyloVersion +"\n\n"
    alignmentReport+=":green:`Method:` ["+ snakemake.config["makeTree"]["method"] + "]_\n\n"
    alignmentReport+="**Command:**\n\n"
    alignmentReport+=":commd:`make_phylogeny.py -i "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned.fasta -o representative_seq_set_noSingletons_aligned_pfiltered.tre "+ snakemake.config["makeTree"]["extra_params"]+ " -t " + snakemake.config["makeTree"]["method"]+"`\n\n"
    alignmentReport+="**Output file:**\n\n"
    alignmentReport+=":green:`- Taxonomy tree:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/aligned/representative_seq_set_noSingletons_aligned.tre\n\n"
    alignmentReport+=makePhyloBenchmark+"\n\n"
#KRONA REPORT
kronaReport = ""
if  snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kronaReport+="Krona report\n----------------\n\n"
    kronaReport+="Krona allows hierarchical data to be explored with zooming, multi-layered pie charts.\n\n"
    kronaReport+=":red:`Tool:` [Krona]_\n\n"
    if snakemake.config["krona"]["otu_table"].casefold() != "singletons":
        kronaReport+="These charts were created using the OTU table **without** singletons\n\n"
    else:
        kronaReport+="These charts were created using the OTU table **including** singletons\n\n"

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
    kronaReport+=" .. _kreport: report_files/krona_report."+snakemake.config["assignTaxonomy"]["tool"]+".html\n\n"

    kronaReport+="Or access the html file at:\n\n"
    kronaReport+=":green:`- Krona html file:` "+snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/otu/taxonomy_"+snakemake.config["assignTaxonomy"]["tool"]+"/krona_report.html\n\n"
    kronaReport+=kronaBenchmark+"\n\n"

###############################################################################
#                         REFERENCES                                     #
################################################################################
#CLUSTER OTUS
if snakemake.config["pickOTU"]["m"] == "uclust":
    variable_refs+= ".. [uclust] Edgar RC. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-2461.\n\n"
elif snakemake.config["pickOTU"]["m"] == "usearch61":
    variable_refs+= ".. [usearch61] Edgar RC. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26(19):2460-2461.\n\n"
elif snakemake.config["pickOTU"]["m"] == "mothur":
    variable_refs+= ".. [mothur] Schloss PD, Wescott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, Lesniewski RA, Oakley BB, Parks DH, Robinson CJ, Sahl JW, Stres B, Thallinger GG, Van Horn DJ, Weber CF. 2009. Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol 75(23):7537-7541.\n\n"
elif snakemake.config["pickOTU"]["m"] == "blast":
    variable_refs+= ".. [blast] Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215(3):403-410\n\n"
elif snakemake.config["pickOTU"]["m"] == "swarm":
    variable_refs+= ".. [swarm] Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2014) Swarm: robust and fast clustering method for amplicon-based studies. PeerJ 2:e593 doi: 10.7717/peerj.593\n\n"
elif snakemake.config["pickOTU"]["m"] == "cdhit":
    variable_refs+= ".. [cdhit] Cd-hit: Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for clustering the next generation sequencing data. Bioinformatics, (2012), 28 (23): 3150-3152. doi: 10.1093/bioinformatics/bts565.\n\n"
#ALIGNMENT
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

report("""
{title}
    .. role:: commd
    .. role:: red
    .. role:: green

**CASCABEL** is designed to run amplicon sequence analysis across single or multiple read libraries. This report consists of the OTU creation and taxonomic assignment for all the combined accepted reads of given samples or libraries, if multiple.

{txtDescription}

Combine Reads
---------------

Merge all the reads of the individual libraries into one single file.

**Command:**

{catCommand}

**Output file:**

:green:`- Merged reads:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta

The total number of reads is: {totalReads}

{combineBenchmark}

{dereplicateReport}

Cluster OTUs
-------------

Assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences based on a user-defined similarity threshold.

:red:`Tool:` [QIIME]_ - pick_otus.py

:red:`Version:` {clusterOtuVersion}

:green:`Method:` [{snakemake.config[pickOTU][m]}]_

:green:`Identity:` {snakemake.config[pickOTU][s]}

**Command:**

:commd:`pick_otus.py -m {snakemake.config[pickOTU][m]} -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/ {snakemake.config[pickOTU][extra_params]} -s {snakemake.config[pickOTU][s]}`

**Output files:**

:green:`- OTU List:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.txt

:green:`- Log file:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.log

The total number of different OTUS is: {totalOtus}

{otuBenchmark}

Pick representatives
-----------------------
Pick a single representative sequence for each OTU.

:red:`Tool:` [QIIME]_ - pick_rep_set.py

:red:`Version:` {pickRepVersion}

:green:`Method:` {snakemake.config[pickRep][m]}

**Command:**

:commd:`pick_rep_set.py -m {snakemake.config[pickRep][m]} -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/seqs_fw_rev_filtered_otus.txt -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/seqs_fw_rev_filtered.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.fasta {snakemake.config[pickRep][extra_params]} --log_fp {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.log`

**Output file:**

:green:`- Fasta file with representative sequences:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/representative_seq_set.fasta

{pikRepBenchmark}

Assign taxonomy
----------------
Given a set of sequences, assign the taxonomy of each sequence.

{assignTaxoStr}

The percentage of successfully assigned OTUs is: {prcAssignedOtus}

**Output file:**

:green:`- OTU taxonomy assignation:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/representative_seq_set_tax_assignments.txt

{assignTaxaBenchmark}

Make OTU table
---------------
Tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU in the last column.

:red:`Tool:` [QIIME]_ - make_otu_table.py

:red:`Version:` {makeOTUVersion}

**Command:**

:commd:`make_otu_table.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/seqs_fw_rev_filtered_otus.txt -t {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/representative_seq_set_tax_assignments.txt {snakemake.config[makeOtu][extra_params]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.biom`

**Output file:**

:green:`- Biom format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.biom

{otuTableBenchmark}

Convert OTU table
------------------
Convert from the BIOM table format to a human readable format.

:red:`Tool:` [BIOM]_

:red:`Version:` {convertBiomVersion}

**Command:**

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.txt {snakemake.config[biom][tableType]} {snakemake.config[biom][headerKey]} {snakemake.config[biom][extra_params]} {snakemake.config[biom][outFormat]}`

**Output file:**

:green:`- TSV format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.txt

{convertOtuBenchmark}

Summarize Taxa
---------------
Summarize information of the representation of taxonomic groups within each sample.

:red:`Tool:` [QIIME]_ - summarize_taxa.py

:red:`Version:` {summTaxaVersion}

**Command:**

:commd:`summarize_taxa.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.biom {snakemake.config[summTaxa][extra_params]} -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/summary/`

**Output file:**

:green:`- Taxonomy summarized counts at different taxonomy levels:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/summary/otuTable_L**N**.txt

Where **N** is the taxonomy level. Default configuration produces levels from 2 to 6.

{summTaxaBenchmark}

Filter OTU table
-----------------
Filter OTUs from an OTU table based on their observed counts or identifier.

:red:`Tool:` [QIIME]_ - filter_otus_from_otu_table.py

:red:`Version:` {filterOTUNoSVersion}

:green:`Minimum observation counts:` {snakemake.config[filterOtu][n]}

**Command:**

:commd:`filter_otus_from_otu_table.py -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable_noSingletons.biom {snakemake.config[filterOtu][extra_params]} -n {snakemake.config[filterOtu][n]}`

**Output file:**

:green:`- Biom table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable_noSingletons.biom

{otuNoSingletonsBenchmark}

Convert Filtered OTU table
---------------------------
Convert the filtered OTU table from the BIOM table format to a human readable format

:red:`Tool:` [BIOM]_

:red:`Version:` {convertBiomVersion}

**Command:**

:commd:`biom convert -i {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable_noSingletons.biom -o {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable_noSingletons.txt {snakemake.config[biom][tableType]} {snakemake.config[biom][headerKey]} {snakemake.config[biom][extra_params]} {snakemake.config[biom][outFormat]}`

**Output file:**

:green:`- TSV format table:` {snakemake.wildcards.PROJECT}/runs/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/otuTable_noSingletons.txt

{otuNoSingletonsBenchmark}

Filter representative sequences
---------------------------------
Remove sequences according to the filtered OTU biom table.

:red:`Tool:` [QIIME]_ - filter_fasta.py

:red:`Version:` {filterFastaVersion}

**Command:**

:commd:`filter_fasta.py -f {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/representative_seq_set.fasta -o {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/representative_seq_set_noSingletons.fasta {snakemake.config[filterFasta][extra_params]} -b {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/otuTable_noSingletons.biom`

**Output file:**

:green:`- Filtered fasta file:` {snakemake.wildcards.PROJECT}/samples/{snakemake.wildcards.run}/otu/taxonomy_{snakemake.config[assignTaxonomy][tool]}/representative_seq_set_noSingletons.fasta

{filterBenchmark}

{alignmentReport}

{kronaReport}

Final counts
-------------

{countTxt}

.. image:: report_files/sequence_numbers_all.png

:red:`Note:`

:green:`- Assigned OTUs percentage` is the amount of successfully assigned OTUs.

:green:`- No singletons percentage` is the percentage of no singletons OTUs in reference to the complete OTU table.

:green:`- Assigned No singletons` is the amount of successfully no singletons assigned OTUs.

References
------------

.. [QIIME] QIIME. Caporaso JG, Kuczynski J, Stombaugh J, Bittinger K, Bushman FD, Costello EK, Fierer N, Gonzalez Pena A, Goodrich JK, Gordon JI, Huttley GA, Kelley ST, Knights D, Koenig JE, Ley RE, Lozupone CA, McDonald D, Muegge BD, Pirrung M, Reeder J, Sevinsky JR, Turnbaugh PJ, Walters WA, Widmann J, Yatsunenko T, Zaneveld J, Knight R. 2010. QIIME allows analysis of high-throughput community sequencing data. Nature Methods 7(5): 335-336.

.. [Cutadapt] Cutadapt v1.15 .Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.Journal, 17(1):10-12, May 2011. http://dx.doi.org/10.14806/ej.17.1.200

.. [vsearch] Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. doi: 10.7717/peerj.2584

.. [Krona] Ondov BD, Bergman NH, and Phillippy AM. Interactive metagenomic visualization in a Web browser. BMC Bioinformatics. 2011 Sep 30; 12(1):385.

.. [BIOM] The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. Daniel McDonald, Jose C. Clemente, Justin Kuczynski, Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob Knight, and J. Gregory Caporaso.GigaScience 2012, 1:7. doi:10.1186/2047-217X-1-7

{variable_refs}


""", snakemake.output[0], metadata="Author: J. Engelmann & A. Abdala ")
