#####
# Dada2- filter and QA reads for Cascabel
# 
# @author: J. Engelmann
# @author: A. ABdala

library(dada2) 
library(Biostrings)

args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = Generate QAplots
#args[3] = trunc FW
#args[4] = trunc RV
#args[5] = maxEE FW
#args[6] = maxEE RV
#args[7] = cpus
#args[8] = filterAndTrim function extra params
#args[9]... = output filter summary
#args[10]... = summary files from libraries
# suppose that we expand and pass the summary.txt from the output
# we should have a path {project}/runs/{run}/{sample}_data/demultiplexed/
#args<-c("/export/lv3/scratch/projects_AA/dada2_CASCABEL/CASCABEL",
#        "T","200","200",3,5,15,",truncQ=2, rm.phix=TRUE",
#       "cascabel_project/runs/Test/asv/filter_summary.out",
#      "cascabel_project/runs/Test/summer_data/demultiplexed/summary.txt", 
#     "cascabel_project/runs/Test/winter_data/demultiplexed/summary.txt")


setwd(args[1])


#Set the different paths for all the supplied libraries
paths = c()
for(i in 10:length(args)) {
  paths <- c(paths,gsub("/summary.txt", '',args[i]))
}

#List files
filesForw <- sort(list.files(paths, pattern="_1.fastq.gz", full.names = TRUE))
filesRev <- sort(list.files(paths, pattern="_2.fastq.gz", full.names = TRUE))

#Get sample names
sample.names <- gsub('_1.fastq.gz', '', basename(filesForw))

#if we want to save the plot we can do the following
asv_dir <- gsub("filter_summary.out"," ", args[9])
dir.create(asv_dir)
if (args[2] == "T" || args[2] == "TRUE" ){
  library(ggplot2)
  plots_fw <- plotQualityProfile(filesForw)
  ggsave(paste0(asv_dir,"fw_QA_plots.pdf"), plots_fw)
  plots_rv <- plotQualityProfile(filesRev)
  ggsave(paste0(asv_dir,"rv_QA_plots.pdf"), plots_fw)
  
  #here we have to solve two things where to store the plots
  # and if we are going to generate one pdf per file, one per library
  # or one per strand...
}

#Create path and file names for filtered samples" 
#filtFs <- file.path(paths, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(paths, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filtFs <- gsub('demultiplexed/','demultiplexed/filtered/',filesForw)
filtRs <- gsub('demultiplexed/','demultiplexed/filtered/',filesRev)

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Within CASCABEL the adapters have been already removed so no extra truncation should be needed
# within  Cacabel we remove barcodes and primers
fw_len <- strtoi(args[3], 10) 
rv_len <- strtoi(args[4], 10)
maxEE_fw <- strtoi(args[5], 10) 
maxEE_rv <- strtoi(args[6],10)
cpus <- strtoi(args[7],10)
extra_params <- args[8]


#The original run looks like this:
#out <- filterAndTrim(filesForw, filtFs, filesRev, filtRs, truncLen=c(fw_len,rv_len),
#                     maxN=0, maxEE=c(maxEE_fw,maxEE_rv), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=5) 

#here we just make it possible to pass extra parameters from the args...

if (!startsWith( trimws(extra_params), ',') && nchar(trimws(extra_params))>1){
  extra_params <- paste(",",extra_params)
}

out <- eval(parse(text=paste("filterAndTrim(filesForw, filtFs, filesRev, filtRs, truncLen=c(fw_len,rv_len),maxN=0, maxEE=c(maxEE_fw,maxEE_rv),compress=TRUE, multithread=cpus",extra_params,")")))   
#The original rownames are the name of the fq files
rownames(out)<-sample.names
colnames(out)<-c("reads.in","filtered")
write.table(out, file = args[9]  , sep = "\t", quote=FALSE, col.names =NA)

