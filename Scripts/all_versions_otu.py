import subprocess
import sys

versions = open(snakemake.output[0],"w")
versions.write("Tool\tVersion\tInformation\n")

########################################################
#                CASCABEL version                      #
########################################################
version = "ND"
try:
  with open('cascabel.version') as f:
    version = f.readline().strip('\n')
except FileNotFoundError:
        print("Cascabel's file version missing: ../cascabel.version\nYou can see README file for Cascabel version.")
        version = "see README file"
versions.write("Cascabel\t"+version+"\tCascabel pipeline\n")


pythonv=str(sys.version_info[0])+"."+str(sys.version_info[1])+"."+str(sys.version_info[2])
versions.write("Python\t"+pythonv+"\tSeveral scripts and snakemake engine.\n")


#--fastq
fqVersion="ND"
if snakemake.config["QC"].lower() ==  "fastqc" or snakemake.config["QC"].lower() ==  "both":
    fqv = subprocess.run([snakemake.config["fastQC"]["command"], '--version'], stdout=subprocess.PIPE)
    fqVersion = fqv.stdout.decode('utf-8').strip() 
    versions.write("FastQC\t"+fqVersion+"\tQuality on fastq files\n")

#--SEQUALI
sequaliVersion = "ND"
if snakemake.config["QC"].lower() ==  "sequali" or snakemake.config["QC"].lower() ==  "both":
    sqv = subprocess.run(['sequali', '--version'], stdout=subprocess.PIPE)
    sequaliVersion = sqv.stdout.decode('utf-8').strip()
    versions.write("Sequali\t"+sequaliVersion+"\tQuality on fastq files\n")
#--Vsearch 
vsearchVersion = "ND"
if snakemake.config["ANALYSIS_TYPE"].lower() ==  "otu":
    vsearchV = subprocess.run(['vsearch', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    vsearchVersion =  vsearchV.stdout.decode('utf-8').split('\n', 1)[0].strip()
    versions.write("Vsearch\t"+vsearchVersion+"\tUsed for: dereplication")
    if snakemake.config["pickOTU"]["m"].lower() ==  "vsearch" :
        versions.write("; clustering")
    if snakemake.config["assignTaxonomy"]["tool"].lower() ==  "vsearch" :
        versions.write("; taxonomy assignation")
    versions.write("\n")

blastnVersion="ND"
if snakemake.config["assignTaxonomy"]["tool"].lower() ==  "blast" :
    blastnV = subprocess.run([snakemake.config["assignTaxonomy"]["blast"]["command"], '-version'], stdout=subprocess.PIPE)
    blastnVersion = blastnV.stdout.decode('utf-8').split('\n', 1)[0].replace('blastn:','').strip()
    versions.write("Blast\t"+blastnVersion+"\tTaxonomy assignation\n")


pearversion = "ND"
try:
    pearv = subprocess.run( [snakemake.config["pear"]["command"]+" -h | grep 'PEAR v'"], stdout=subprocess.PIPE, shell=True)
    pearversion =  pearv.stdout.decode('utf-8').strip() 
except Exception as e:
    pearversion = "ND"
versions.write("Pear\t"+pearversion+"\tAssembly (pair) reads\n")

#--cutadapt
cutVersion = "ND"
if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg" or snakemake.config["demultiplexing"]["demultiplex"].lower() == "t":
    cutv = subprocess.run(['cutadapt', '--version'], stdout=subprocess.PIPE)
    cutVersion =  cutv.stdout.decode('utf-8').strip()
    versions.write("Cutadapt\t"+cutVersion+"\t")
    if snakemake.config["primers"]["remove"].casefold() == "metadata" or snakemake.config["primers"]["remove"].casefold() == "cfg":
        versions.write("Remove primer;")
    if snakemake.config["demultiplexing"]["demultiplex"].lower() == "t":
        versions.write(" Demultiplex")
    versions.write("\n") 

biomVersion="ND"
biomV = subprocess.run([snakemake.config["biom"]["command"], '--version'], stdout=subprocess.PIPE)
biomVersion =  biomV.stdout.decode('utf-8').strip() 
versions.write("Biom\t"+biomVersion+"\tConvert tables to biom format.\n")

if  snakemake.config["pickOTU"]["m"] == "swarm":
    swarmVersion="ND"
    swarmV = subprocess.run(['swarm', '--version'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    swarmVersion = "**" + vsearchV.stdout.decode('utf-8').split('\n', 1)[0].strip() + "**"
    versions.write("Swarm\t"+swarmVersion+"\tOTU Clustering.\n")
    

if  snakemake.config["krona"]["report"].casefold() == "t" or snakemake.config["krona"]["report"].casefold() == "true":
    kVersion="ND"
    try:
        kv=subprocess.run([snakemake.config["krona"]["command"]], stdout=subprocess.PIPE)
        kVersion=kv.stdout.decode('utf-8').split('\n')[1].strip().replace("_","").replace("/","").replace("\\","").split("-")[0].strip()
    except Exception as e:
        kVersion="ND"
    versions.write("Krona Tools\t"+kVersion+"\tKrona interactive report.\n")

