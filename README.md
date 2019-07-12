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
+ OTU tables
+ Taxonomy OTU assignation
+ Taxonomy summary
+ Representative sequence alignment
+ Phylogenetic tree
+ CASCABEL Report

## Install the pipeline

You can clone the project

<pre><code class="text">
git clone https://github.com/AlejandroAb/CASCABEL.git
</code></pre>

Download it from this repository

<pre><code class="text">
https://github.com/AlejandroAb/CASCABEL
</code></pre>


## Load Snakemake environment

In order to have snakemake running into your shell, you need to be located at your $HOME directory (/export/data/username)
and execute the following command:

<pre><code class="text">
source  .zshrc.conda344
</code></pre>


To return to your original shell environment just type:

<pre><code class="php">
source  .zshrc
</code></pre>

## Configure pipeline

In order to use any of these implementations, just copy the complete subdirectory into your personal directories and read further instructions at the manual.pdf located at /export/data
document for the specific analysis you are planning to perform.

For a complete guide and complete documentation please visit the official project [wiki](http://redmine.nioz.nl/projects/pipeline-for-amplicon-analysis/wiki/Wiki) 
