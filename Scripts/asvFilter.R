#####
# Dada2- filter and QA reads for Cascabel
# 
# @author: J. Engelmann
# @author: A. Abdala

library(dada2) 
library(Biostrings)

args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = Generate QAplots (T or F)
#args[3] = trunc FW 
#args[4] = trunc RV
#args[5] = maxEE FW
#args[6] = maxEE RV
#args[7] = cpus
#args[8] = filterAndTrim function extra params
#args[9] = interactive behavior
#args[10] = output filter summary
#args[11] = remove primers: config["demultiplexing"]["primers"]["remove"]
#args[12]... = summary files from libraries


setwd(args[1])

#Set the different paths for all the supplied libraries
paths = c()
if (args[11] != "F" && args[11] != "FALSE" ){
  for(i in 12:length(args)) {
    paths <- c(paths,gsub("/summary.txt", '/primer_removed',args[i]))
  }
}else{
 for(i in 12:length(args)) {
    paths <- c(paths,gsub("/summary.txt", '',args[i]))
  }
}

#List files
filesForw <- sort(list.files(paths, pattern="_1.fastq.gz$", full.names = TRUE))
filesRev <- sort(list.files(paths, pattern="_2.fastq.gz$", full.names = TRUE))

#Get sample names
sample.names <- gsub('_1.fastq.gz', '', basename(filesForw))

print(sample.names)

#if we want to save the plot we can do the following
asv_dir <- gsub("filter_summary.out","", args[10])
dir.create(asv_dir, showWarnings = TRUE)
if (args[2] == "T" || args[2] == "TRUE" ){
  library(ggplot2)
  plots<-list()
  for(i in 1:length(filesForw)){
    p<-plotQualityProfile(c(filesForw[i],filesRev[i]))
    plots<-append(plots,list(p))
  }
  pdf(paste0(asv_dir,"QA_plots.pdf"))
  print(plots)
  dev.off()
}

readintegerConsole <- function(x){
  cat(x);
  n <- readLines("stdin",n=1);
  if(!grepl("^[0-9]+$",n))
  {
    return(readintegerConsole(x))
  }
  return(as.integer(n))
}

createMenuConsole <- function()
{ 
  cat("\n--------------------------------------\n")
  cat("------     READ TRUNCATION    --------\n")
  cat("The next step will truncate the reads:\n")
  if (args[2] == "T" || args[2] == "TRUE" ){
   cat(paste("We advise having a look at the QA plots generated at: ",asv_dir))
  }else{
    cat("Based on your configuration settings, no new QC plots were generated.\n")
    cat("In order to select the best values, we advise you to have a look at your\n")
    cat("FastQC report but be aware that no barcodes or primers were removed when\n")
    cat("they were generated, so probably you have shorter reads at this stage.\n") 
    cat("If you want to generate dada2's QA plots, stop the workflow with option 3\n")
    cat("then update your settings with \"dada2_filter: generateQAplots = T\" and \n")
    cat("re-run the pipeline.\n\n") 
  }
  cat("\nYour current values are:\n")
  cat(paste(" * Forward: ",args[3],"\n"))
  cat(paste(" * Reverse: ",args[4],"\n"))
  cat("\nA value of \"0\" means no truncation.\n")
  cat("Reads shorter than the selected value are discarded!\n")
  cat("What would you like to do?\n\n")
  
  cat(" 1. Use values from configuration file.\n")
  cat(" 2. Specify new values!\n")
  cat(" 3. Interrupt workflow.\n")
  cat("Enter your option:\n")
  n <- readLines("stdin",n=1);
  if(!grepl("^[0-9]+$",n))
  {
    return(createMenuConsole())
  }else if(n == 1){
    fw_len <<- strtoi(args[3], 10) 
    rv_len <<- strtoi(args[4], 10)
  }else if(n == 2){
    
    fw_len <<- readintegerConsole("Enter FW reads truncation value:")
    rv_len <<- readintegerConsole("Enter RV reads truncation value:")
  
  }else if(n == 3){
    exit(1)
  }else {
    return(createMenuConsole())
  }
}

#Create path and file names for filtered samples" 
#filtFs <- file.path(paths, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(paths, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

if (args[11] != "F" && args[11] != "FALSE" ){
  filtFs <- gsub('primer_removed/','filtered/',filesForw)
  filtRs <- gsub('primer_removed/','filtered/',filesRev)

}else{
  filtFs <- gsub('demultiplexed/','demultiplexed/filtered/',filesForw)
  filtRs <- gsub('demultiplexed/','demultiplexed/filtered/',filesRev)
}

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Within CASCABEL the adapters have been already removed so no extra truncation should be needed
# within  Cacabel we remove barcodes and primers
if (args[9] == "T" || args[9] == "TRUE" ){
  createMenuConsole()
}else{
  fw_len <- strtoi(args[3], 10) 
  rv_len <- strtoi(args[4], 10)
  
}
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
write.table(out, file = args[10]  , sep = "\t", quote=FALSE, col.names =NA)

