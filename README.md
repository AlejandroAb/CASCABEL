## CASCABEL

**CASCABEL** is a pipeline designed to run amplicon sequence analysis across single or multiple read libraries.
The objective of this pipeline is to create different output files which allow the user to explore data in a simple and
meaningful way, as well as facilitate downstream analysis, based on the generated output files.

CASCABEL was designed for short read high-throughput sequence data. It covers quality control on the fastq files, assembling paired-end reads to fragments (it can also handle single end data), splitting the libraries into samples (optional), OTU picking and taxonomy assignment. Besides other output files, it will return an OTU table.

Our pipeline is implemented with Snakemake as workflow management engine and allows customizing the analyses by offering several choices for most of the steps. The pipeline can make use of multiple computing nodes and scales from personal computers to computing servers. The analyses and results are fully reproducible and documented in an html and optional pdf report.

## Installation

The easiest and recommended way to do install **Cascabel**  is via **Conda**. The fastest way to obtain Conda is to install Miniconda, a mini version of Anaconda that includes only conda and its dependencies. 

### Miniconda

In order to install conda or miniconda please see the following [tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) (recommended) or, if you are working with a Linux OS, you can try the following:

Download the installer:
<pre><code class="text">
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
</code></pre>

Execute the installation script and follow the instructions.
<pre><code class="text">
bash Miniconda3-latest-Linux-x86_64.sh
</code></pre>

Unfortunately Cascabel have many dependencies and latest Conda releases find conflicts among them, however with conda v 4.6.14 we noticed that the installation can run smoothly. In order to do so, we need to downgrade the conda version with the following command: 

<pre><code class="text">
conda install conda=4.6.14
</code></pre>

### Download CASCABEL

Once that you have **conda** installed we are ready to clone or download the project.

You can clone the project:

<pre><code class="text">
git clone https://github.com/AlejandroAb/CASCABEL.git
</code></pre>

Or download it from this repository:

<pre><code class="text">
wget https://github.com/AlejandroAb/CASCABEL/archive/master.zip
</code></pre>

After downloading or cloning the repository, cd to the "CASCABEL" directory and there execute the following command in order to create CASCABEL's environment:

<pre><code class="text">
conda env create --name cascabel --file environment.yaml
</code></pre>

### Snakemake

Now that you have cascabel's environment created, you can install Snakemake following this [on line help](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or execute the following command:

<pre><code class="text">
conda install -c bioconda -c conda-forge snakemake
</code></pre>

### Matplotlib
 
All the dependencies required by CASCABEL except by Snakemake and thus Python are loaded in one [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html). In this sense, CASCABEL uses [matplotlib](https://matplotlib.org/) for generating some charts, therefore it is needed to have this library installed prior to load the environment. The recommended way to do this is following the [installation guide](https://matplotlib.org/users/installing.html) or you can also try with:
<pre><code class="text">
pip install matplotlib  --user
</code></pre>

*Consider to use the flag *--user* as above, if you are doing a local installation or if you don't have *sudo* rights

### Activate environment

After installing Snakemake and Matplotlib we can activate our new environment. 

<pre><code class="text">
conda activate cascabel
</code></pre>

_After activating the environment it is possible that Snakemake is not more in your PATH, in such case just export Snakemake's bin directory. i.e:_

<pre><code class="text">
export PATH=$PATH:/path/to/miniconda3/bin
</code></pre>

### DADA2

You only need to follow this one more step if you are planning to run **Cascabel** with the *asv workflow*.

There are some issues reported while installing dada2 within conda, therefore we need to perform one more final step in order to install **dada2**

Enter to R shell (just type `R`) and execute the following command:  

<pre><code class="text">
BiocManager::install("dada2", version = "3.10")
</code></pre>

*Please notice that BiocManager should be already installed, so you just need to execute previous command. You can also find more information at [dada2's installation guide.](https://benjjneb.github.io/dada2/dada-installation.html)

### Singularity

We are aware that this is not the easiest installation, therefore we are working on a singularity container, same that we hope to have available soon.

Thanks for your understanding! 

## Getting started

**Required input files:**
                         
+ Forward raw reads (fastq or fastq.gz)
+ Reverse raw reads (fastq or fastq.gz)
+ File with barcode information (only for demultiplexing: [format](http://qiime.org/documentation/file_formats.html#metadata-mapping-files))

  
**Main expected output files for downstream analysis**

+ Demultiplexed and trimmed reads
+ OTU table
+ Representative sequences fasta file
+ Taxonomy OTU assignation
+ Taxonomy summary
+ Representative sequence alignment
+ Phylogenetic tree
+ CASCABEL Report

**Run CASCABEL**

All the parameters and behavior of the workflow is specified through the [configuration file](../../wiki#23-the-configuration-file-configyaml), therefore the easiest way to have the pipeline running is to filling up some required parameters on such file.

<pre><code class="yaml">
#------------------------------------------------------------------------------#
#                             Project Name                                     #
#------------------------------------------------------------------------------#
# The name of the project for which the pipeline will be executed. This should #
# be the same name used as the first parameter on init_sample.sh script (if    #
# used for multiple libraries                                                 #
#------------------------------------------------------------------------------#
PROJECT: "My_CASCABEL_Project"

#------------------------------------------------------------------------------#
#                            LIBRARIES/SAMPLES                                 #
#------------------------------------------------------------------------------#
# SAMPLES/Libraries you will like to include on the analysis                   #
# Same library names used  with init_sample.sh script (if used)                #
# Include each sample / library name surrounded by quotes, and comma separated #
# i.e LIBRARY["SAMP1","SAMP2",..."SAMPN"]                                      #
#------------------------------------------------------------------------------#
LIBRARY: ["EXP1"]
#------------------------------------------------------------------------------#
#                             INPUT FILES                                      #
#------------------------------------------------------------------------------#
# Here you have to enter the FULL PATH for both: the raw reads and the metadata# 
# file (barcode mapping file). The metadata file is only needed if you want to #
# perform demultiplexing.                                                      #
# - fw_reads:  Full path to the raw reads in forward direction (R1)            #
# - rw_reads:  Full path to the raw reads in reverse direction (R2)            #
# - metadata:  Full path to the metadata file with barcodes for each sample    #
#              to perform library demultiplexing                               #
#------------------------------------------------------------------------------#
fw_reads: "/full/path/to/forward.reads.fq"
rv_reads: "/full/path/to/reverse.reads.fq"
metadata: "/full/path/to/metadata.barcodes.txt"

#------------------------------------------------------------------------------#
#                               RUN                                            #
#------------------------------------------------------------------------------#
# Name of the RUN - Only use alphanumeric characters and don't use spaces.     #
# This parameter helps the user to execute different runs (pipeline executions)#
# with the same input data but with different parameters (ideally).            #
# The RUN parameter can be set here or remain empty, in the latter case, the   #
# user must assign this value via the command line.                            #
# i.e:  --config RUN=run_name                                                  #
#------------------------------------------------------------------------------#
RUN: "My_First_run"

#------------------------------------------------------------------------------#
#                                 ANALYSIS TYPE                                #
# rules:                                                                       #
#------------------------------------------------------------------------------#
# Cascabel supports two main types of analysis:                                #
#  1) Analysis based on traditional OTUs (Operational Taxonomic Units) which   #
#     are mainly generated by clustering sequences based on a sheared          #
#     similarity threshold.                                                    #
#  2) Analysis based on ASV (Amplicon sequence variant). This kind of analysis #
#     deal also with the errors on the sequence reads such that true sequence  #
#     variants can be resolved, down to the level of single-nucleotide         #
#     differences.                                                             #
#                                                                              #
#-----------------------------       PARAMS       -----------------------------#
#                                                                              #
# - ANALYSIS_TYPE    "OTU" or "ASV". Defines the type analysis                 #
#------------------------------------------------------------------------------#
ANALYSIS_TYPE: "OTU"


</code></pre>

_If you need to run CASCABEL for multiple libraries please follow this [instructions](../../wiki#22-initialize-structure-for-multiple-libraries)_

As you can see on the previous fragment of the configuration file (config.yaml), the required parameters for CASCABEL to start are: *PROJECT*, *LIBRARY*, *RUN*,  *fw_reads*, *rv_reads* and *metadata*. After entering these parameters, take some minutes and go through the rest of the config file and overwrite settings according to your needs. Most values are already pre-configured. The config file explains itself by using meaningful headers before each rule, explaining the aim of such rule and the different parameters the user can use. It is very important to keep the indentation of the file (don’t change the tabs
and spaces), as well as the name of the parameters. Once that you have valid values for these entries, you are ready to run the pipeline (before start CASCABEL always is a good practice to make a ["dry run"](../../wiki#31-dry-run)):

Also, please notice the **ANALYSIS TYPE** section. Cascabel, supports two main type of analysis, OTUs (Operational Taxonomic Units) and ASVs (Amplicon Sequence Variants), here you can select the target workflow that Cascabel will execute. For more information pleasee refer to the [Analysis type section](../../wiki#233-analysis-type) 
 
<pre><code class="yaml">
snakemake --configfile config.yaml 
</code></pre>

Optionally you can specify the same parameters* via --config flag, rather than within the config.yaml file:

<pre><code class="text">
 snakemake --configfile config.yaml --config PROJECT="My_CASCABEL_Project"  RUN="My_First_run" fw_reads="//full/path/to/forward.reads.fq" rv_reads="/full/path/to/reverse.reads.fq" metadata="full/path/to/metadata.barcodes.txt"
</code></pre>


*Except for the LIBRARY, as this is declared as an array, therefore it must be filled up within the configuration file

### Configure pipeline

For a complete guide on how to setup and use CASCABEL please visit the official project [wiki](https://github.com/AlejandroAb/CASCABEL/wiki)

### Configuration files

We supply some "pre-filled" configuration files for the main possible configurations like for double and single barcoded paired end reads for OTU and ASV analysis. We strongly advise to make informed choices about parameter settings matching the individual needs of the experiment and data set.

### Test data

In order to test the pipeline we also sugest to try running it with [CASCABEL's test data](https://github.com/AlejandroAb/CASCABEL-Test)

### Citing

For the moment Cascabel's prepint is available for citing at [www.biorxiv.org](https://www.biorxiv.org/content/early/2019/10/17/809384)

Cascabel: a flexible, scalable and easy-to-use amplicon sequence data analysis pipeline
Alejandro Abdala Asbun, Marc A Besseling, Sergio Balzano, Judith van Bleijswijk, Harry Witte, Laura Villanueva, Julia C Engelmann
bioRxiv 809384; doi: https://doi.org/10.1101/809384 
