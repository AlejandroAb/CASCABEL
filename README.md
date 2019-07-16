## CASCABEL

**CASCABEL** is a pipeline designed to run amplicon sequence analysis across single or multiple read libraries.
The objective of this pipeline is to create different output files which allow the user to explore data in a simple and
meaningful way, as well as facilitate downstream analysis, based on the generated output files.


## Getting started

This pipeline is intended to provide a simple and comprehensive work-flow for amplicon sequence data analysis. The pipeline creates different output files which allow the user to explore the data and results in a simple way, as well as facilitate downstream analysis based on the generated output files.

CASCABEL was designed for short read high-throughput sequence data. It covers quality control on the fastq files, assembling paired-end reads to fragments (it can also handle single end data), splitting the libraries into samples (optional), OTU picking and taxonomy assignment. Besides other output files, it will return an OTU table.
The pipeline can be customized to work on different marker genes, reference data bases and applying a range of different methods for e.g. picking OTUs. Many steps are based on QIIME scripts and the full functionality of these scripts is available to the user.

                      
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

## Install the pipeline

The pipeline relays on Snakemake as the workflow-engine, therefore, it is necessary to install this tool in order to have the pipeline running. The easiest and recommended way to install Snakameke is via Conda. The fastest way to obtain Conda is to install Miniconda, a mini version of Anaconda that includes only conda and its dependencies. 

In order to install conda or miniconda please see the following [tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html)

Now that you have conda installed, you can install Snakemake following this on-line help or just with the following command:

<pre><code class="text">
conda install -c bioconda -c conda-forge snakemake
</code></pre>

Once that you have Snakemake installed we are ready to clone or download the project.

You can clone the project

<pre><code class="text">
git clone https://github.com/AlejandroAb/CASCABEL.git
</code></pre>

Download it from this repository

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


## Configure pipeline

For a complete guide and complete documentation please visit the official project [wiki](http://redmine.nioz.nl/projects/pipeline-for-amplicon-analysis/wiki/Wiki) 
