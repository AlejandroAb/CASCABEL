## CASCABEL

**CASCABEL** is a pipeline designed to run amplicon sequence analysis across single or multiple read libraries.
The objective of this pipeline is to create different output files which allow the user to explore data in a simple and
meaningful way, as well as facilitate downstream analysis, based on the generated output files.

CASCABEL was designed for short read high-throughput sequence data. It covers quality control on the fastq files, assembling paired-end reads to fragments (it can also handle single end data), splitting the libraries into samples (optional), OTU picking and taxonomy assignment. Besides other output files, it will return an OTU table.

Our pipeline is implemented with Snakeme as workflow management engine and allows customizing the analyses by offering several choices for most of the steps. The pipeline can make use of multiple computing nodes and scales from personal computers to computing servers. The analyses and results are fully reproducible and documented in an html and optional pdf report.

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

Download it from this repository:

<pre><code class="text">
https://github.com/AlejandroAb/CASCABEL
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
+ Taxonomy OTU assignation
+ Taxonomy summary
+ Representative sequence alignment
+ Phylogenetic tree
+ CASCABEL Report

**Run CASCABEL**

All the parameters and behavior of the workflow is specified through the configuration file, therefore the easiest way to have the pipeline running is to filling up some required variables on such file: 

<pre><code class="yaml">
#------------------------------------------------------------------------------#
#                             Project Name                                     #
#------------------------------------------------------------------------------#
# The name of the project for which the pipe line will be executed. This should#
# be the same name used as the first parameter on init_sample.sh script        #
#------------------------------------------------------------------------------#
PROJECT: "My_CASCABEL_Project"

#------------------------------------------------------------------------------#
#                            LIBRAIRES/SAMPLES                                 #
#------------------------------------------------------------------------------#
# SAMPLES/Libraries you will like to include on the analysis                   #
# Same library names used  with init_sample.sh script                          #
# Include each sample name surrounded by quotes, and comma separated           #
# i.e LIBRARY["SAMP1","SAMP2",..."SAMPN"]                                      #
#------------------------------------------------------------------------------#
LIBRARY: ["EXP1"]
#------------------------------------------------------------------------------#
#                             INPUT FILES                                      #
#------------------------------------------------------------------------------#
fw_reads: "/full/path/to/forward.reads.fq"
rv_reads: "/full/path/to/reverse.reads.fq"
metadata: "/full/path/to/metadata.barcodes.txt"

#------------------------------------------------------------------------------#
#                               RUN                                            #
#------------------------------------------------------------------------------#
# Name of the RUN - Only use alphanumeric characters and don't use spaces.     #
# This argument helps the user to execute different runs (pipeline execution)  #
# with the same input data but with different parameters (ideally).            #
# The RUN variable can be set here or remain empty, in the latter case, the    #
# user must assign this value via the command line --config RUN=User_run_name  #
#------------------------------------------------------------------------------#
RUN: "My_First_run"

</code></pre>
As you can see on the previous config.yaml file the variables that are required for CASCABEL to start are:
PROJECT, LIBRARY, RUN,  fw_reads, rv_reads and metadata. Once that you have valid values for these entries, you are ready to run the pipeline:

<pre><code class="yaml">
snakemake --configfile config.yaml 
</code></pre>

Optionally you can also enter the same variables* via --config flag:

<pre><code class="text">
 snakemake --configfile config.yaml --config PROJECT="My_CASCABEL_Project"  RUN="My_First_run" fw_reads="//full/path/to/forward.reads.fq" rv_reads="/full/path/to/reverse.reads.fq" metadata="full/path/to/metadata.barcodes.txt"
 </code></pre>


* Except for the LIBRARY, as this is declared as an array, therefore it must be filled up within the config.yaml file

### Configure pipeline

For a complete guide on how to setup and use CASCABEL please visit the official project [wiki](http://redmine.nioz.nl/projects/pipeline-for-amplicon-analysis/wiki/Wiki) 
