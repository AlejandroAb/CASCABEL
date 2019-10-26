## CASCABEL

**CASCABEL** is a pipeline designed to run amplicon sequence analysis across single or multiple read libraries.
The objective of this pipeline is to create different output files which allow the user to explore data in a simple and
meaningful way, as well as facilitate downstream analysis, based on the generated output files.

CASCABEL was designed for short read high-throughput sequence data. It covers quality control on the fastq files, assembling paired-end reads to fragments (it can also handle single end data), splitting the libraries into samples (optional), OTU picking and taxonomy assignment. Besides other output files, it will return an OTU table.

Our pipeline is implemented with Snakemake as workflow management engine and allows customizing the analyses by offering several choices for most of the steps. The pipeline can make use of multiple computing nodes and scales from personal computers to computing servers. The analyses and results are fully reproducible and documented in an html and optional pdf report.

## Installation

The first requirement to have the pipeline running, is to install Snakemake. The easiest and recommended way to do this is via Conda. The fastest way to obtain Conda is to install Miniconda, a mini version of Anaconda that includes only conda and its dependencies. 

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

### Snakemake


Now that you have conda installed, you can install Snakemake following this [on line help](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) or execute the following command:

<pre><code class="text">
conda install -c bioconda -c conda-forge snakemake
</code></pre>

### Matplotlib
 
All the dependencies required by CASCABEL except by Snakemake and thus Python are loaded in one [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html). In this sense, CASCABEL uses [matplotlib](https://matplotlib.org/) for generating some charts, therefore it is needed to have this library installed prior to load the environment. The recommended way to do this is following the [installation guide](https://matplotlib.org/users/installing.html) or you can also try with:
<pre><code class="text">
pip install matplotlib
</code></pre>

### Download CASCABEL

Once that you have Snakemake installed we are ready to clone or download the project.

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

Activate the environment:

<pre><code class="text">
conda activate cascabel
</code></pre>

*After activating the environment is possible that Snakemake is not more in your PATH, in such case just export Snakemake's bin directory. i.e:*

<pre><code class="text">
export PATH=$PATH:/path/to/miniconda3/bin
</code></pre>

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

</code></pre>

_If you need to run CASCABEL for multiple libraries please follow this instructions(../../wiki#22-initialize-structure-for-multiple-libraries)_

As you can see on the previous fragment of the configuration file (config.yaml), the required parameters for CASCABEL to start are: *PROJECT*, *LIBRARY*, *RUN*,  *fw_reads*, *rv_reads* and *metadata*. After entering these parameters, take some minutes and go through the rest of the config file and overwrite settings according to your needs. Most values are already pre-configured. The config file explains itself by using meaningful headers before each rule, explaining the aim of such rule and the different parameters the user can use. It is very important to keep the indentation of the file (don’t change the tabs
and spaces), as well as the name of the parameters. Once that you have valid values for these entries, you are ready to run the pipeline (before start CASCABEL always is a good practice to make a ["dry run"](../../wiki#31-dry-run)):

<pre><code class="yaml">
snakemake --configfile config.yaml 
</code></pre>

Optionally you can specify the same parameters* via --config flag, rather than within the config.yaml file:

<pre><code class="text">
 snakemake --configfile config.yaml --config PROJECT="My_CASCABEL_Project"  RUN="My_First_run" fw_reads="//full/path/to/forward.reads.fq" rv_reads="/full/path/to/reverse.reads.fq" metadata="full/path/to/metadata.barcodes.txt"
 </code></pre>


*Except for the LIBRARY, as this is declared as an array, therefore it must be filled up within the config.yaml file

### Configure pipeline

For a complete guide on how to setup and use CASCABEL please visit the official project [wiki](https://github.com/AlejandroAb/CASCABEL/wiki)

### Test data

In order to test the pipeline we also sugest to try running it with [CASCABEL's test data](https://github.com/AlejandroAb/CASCABEL-Test)

### Citing

For the moment Cascabel's prepint is available for citing at [www.biorxiv.org](https://www.biorxiv.org/content/early/2019/10/17/809384)

Cascabel: a flexible, scalable and easy-to-use amplicon sequence data analysis pipeline
Alejandro Abdala Asbun, Marc A Besseling, Sergio Balzano, Judith van Bleijswijk, Harry Witte, Laura Villanueva, Julia C Engelmann
bioRxiv 809384; doi: https://doi.org/10.1101/809384 
