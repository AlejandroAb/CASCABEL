
#####
# Dada2- learn errors and run dada2
#
# @author: J. Engelmann
# @author: A. ABdala

library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = T)
#params <- snakemake@params
#args <- unlist((strsplit(unlist(params), split=" ")))

#args[1] = path for seting WD to use current dir use $PWD
#args[2] = dada2 pool option
#args[3] = cpus
#args[4] = plotErr
#args[5] = extra parameters for dada2 function
#args[6] = path for outputfiles
#args[7] = shorts
#args[8] = longs
#args[9] = ofsets
#args[10]... = chimera
#args[11]... = taxa_db_path
#args[12]... = species_db_path
#args[13]... = add_species
#args[14]... = extra_paramas taxonomy
#args[15]... = minOverlap merge
#args[16]... = maxMismatch merge
#args[17]... = summary files from libraries

# suppose that we expand and pass the summary.txt from the output
# we should have a path {project}/runs/{run}/{sample}_data/demultiplexed/
#args<-c("/export/lv3/scratch/projects_AA/dada2_CASCABEL/CASCABEL",
 #       "pseudo","5","T","to_remove",
  #      "cascabel_project/runs/Test/asv/","250","255","10","T",
   #     "/export/data01/databases/silva/r132/dada2/silva_nr_v132_train_set.fa.gz",
    #   "/export/data01/databases/silva/r132/dada2/silva_species_assignment_v132.fa.gz",
     #   "T",",minBoot=50",12,0,
      #  "cascabel_project/runs/Test/summer_data/demultiplexed/summary.txt", 
       # "cascabel_project/runs/Test/winter_data/demultiplexed/summary.txt","cascabel_project/runs/Test")
#print(args)
setwd(args[1])
#print(getwd())
#Set the different paths for all the supplied libraries
#print(args[1])
#paths = c()
#with args
paths <-NULL
for(i in 17:(length(args)-1)) {
  paths <- c(paths,gsub("/summary.txt", '',args[i]))
}
print(paths)
#with params
#inputs <- unlist((strsplit(unlist(snakemake@input[[2]]), split=" ")))
#print(inputs)
#for(i in 1:length(inputs)) {
#    paths <- c(paths,gsub("/summary.txt", '',inputs[i]))
#}
#print(paths)
#List files
filesForw <- sort(list.files(paths, pattern="_1.fastq.gz", full.names = TRUE))
filesRev <- sort(list.files(paths, pattern="_2.fastq.gz", full.names = TRUE))
#Get sample names
sample.names <- gsub('_1.fastq.gz', '', basename(filesForw))
#Create path and file names for filtered samples"
#filtFs <- file.path(paths, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
#filtRs <- file.path(paths, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filtFs <- gsub('demultiplexed/','demultiplexed/filtered/',filesForw)
filtRs <- gsub('demultiplexed/','demultiplexed/filtered/',filesRev)
#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if (args[2] == "pseudo"){
  pool = "pseudo"
}else{
  pool <- eval(parse(text=args[1]))
}
cpus <- strtoi(args[3],10)
extra_params <- args[5]
# learn error rates
errF <- learnErrors(filtFs, multithread=cpus)
#errF <- eval(parse(text=paste("learnErrors(filtFs, multithread=cpus,",extra_params,")")))
errR <- learnErrors(filtRs, multithread=cpus)
#errR <- eval(parse(text=paste("learnErrors(filtRs, multithread=cpus)")))
#plotErrors(errF, nominalQ=TRUE)
if (args[4] == "T" || args[4] == "TRUE" ){
  library(ggplot2)
  plots_fw <- plotErrors(errF, nominalQ=TRUE)
  ggsave(paste(args[6],"fr_err.pdf",sep=""), plots_fw)
  plots_rv <- plotErrors(errR, nominalQ=TRUE)
  ggsave(paste(args[6],"rv_err.pdf", sep=""), plots_rv)
  #here we have to solve two things where to store the plots
  # and if we are going to generate one pdf per file, one per library
  # or one per strand...
}



dadaFs <- dada(filtFs, err=errF, multithread=cpus, pool=pool)
dadaRs <- dada(filtRs, err=errR, multithread=cpus, pool=pool)

#if (!startsWith( trimws(extra_params), ',') && nchar(trimws(extra_params))>1){
#  extra_params <- paste(",",extra_params)
#}
#dadaFs <- eval(parse(text=paste("dada(filtFs, err=errF, multithread=cpus, pool=pool, ",extra_params,")")))   
#dadaRs <- eval(parse(text=paste("dada(filtRs, err=errF, multithread=cpus, pool=pool, ",extra_params,")")))   

#merge
minOv <- as.integer(args[15])
maxMism <- as.integer(args[16])
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = minOv, maxMismatch = maxMism)

# ASV table 
seqtab <- makeSequenceTable(mergers)
#rownames= samples #colnames= sequences

#this show sequence length distributions
seq_hist <- table(nchar(getSequences(seqtab)))
fname_seqh <- paste(args[6],"seq_hist.txt",sep="")
write.table(seq_hist, file = fname_seqh  , sep = "\t", quote=FALSE, col.names = FALSE)

#fname_asv_obj <- paste(args[6],"asv.rds")
# Save an object to a file
#saveRDS(object, file = fname_asv_obj)

shorts <- as.integer(args[7])
longs <- as.integer(args[8])
offset <- as.integer(args[9])
m <- mean(nchar(getSequences(seqtab)))

readinteger <- function(x)
{ 
  n <- readline(prompt=x)
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger(x))
  }
  
  return(as.integer(n))
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
  cat("Please enter the option which fits better for your data:\n")
  cat(paste("1. Use values from configuration file: length > ", shorts," and length < ", longs,"\n"))
  cat(paste("2. Use values from median + /- the offset: length > ", (m-offset)," and length < ", (m+offset),"\n"))
  cat("3. Specify new values!\n")
  cat("4. Print sequence length histogram.\n")
  cat("5. Interrupt workflow.\n")
  cat("Enter your option:\n")
  n <- readLines("stdin",n=1);
  if(!grepl("^[0-9]+$",n))
  {
    return(createMenuConsole())
  }else if(n == 1){
    shorts <<- as.integer(args[7])
    longs <<- as.integer(args[8])
  }else if(n == 2){
    shorts <<- as.integer(m-offset)
    longs  <<- as.integer(m +offset)
  }else if(n == 3){
    
    shorts <<- readintegerConsole("Enter the shortes length allowed:")
    longs <<- readintegerConsole("Enter the longest length allowed:")
  }
  else if(n == 4){
    print(seq_hist)
    return(createMenuConsole())
  }else if(n == 5){
    exit(1)
  }else {
    return(createMenuConsole())
  }
}
createMenu <- function()
{ 
  print("Please enter the option which fits better for your data:")
  print(paste("1. Use values from configuration file: length > ", shorts," and length < ", longs))
  print(paste("2. Use values from median + /- the offset: length > ", (m-offset)," and length < ", (m+offset)))
  print("3. Specify new values!")
  print("4. Print sequence length histogram.")
  print("5. Interrupt workflow.")
  n <- readline(prompt="Enter your option:") 
  if(!grepl("^[0-9]+$",n))
  {
    return(createMenu())
  }else if(n == 1){
    shorts <<- as.integer(args[7])
    longs <<- as.integer(args[8])
  }else if(n == 2){
    shorts <<- as.integer((m-offset))
    longs  <<- as.integer((m +offset))
  }else if(n == 3){
    
    shorts <<- readinteger("Enter the shortes length allowed:")
    longs <<- readinteger("Enter the longest length allowed:")
  }
  else if(n == 4){
    print(seq_hist)
    return(createMenu())
  }else if(n == 5){
    exit(1)
  }else {
    return(createMenu())
  }
}

createMenuConsole()
print(shorts)
print(longs)
write.table(c(shorts,longs),paste0(args[6],"shorts_longs.log"), quote=F,sep="\t", row.names = F, col.names = F)
#This is analogous to “cutting a band” in-silico to get amplicons of the targeted length
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% shorts:longs] 

if (args[10] == "T" || args[10] == "TRUE" ){
   
# rm chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=cpus, verbose=TRUE)  
  #dim(seqtab.nochim)
  #sum(seqtab.nochim)/sum(seqtab)
  new_names <- c(paste("asv.",1:length(colnames(seqtab.nochim)),sep=""))
  seqs <- getSequences(seqtab.nochim)
  #estos son equivalentes colnames(seqtab.nochim) == seqs
  #names(seqs) <- seqs
  names(seqs) <- new_names
  original_names<-colnames(seqtab.nochim)
  colnames(seqtab.nochim)<-new_names
  # also write to file
  write.table(seqtab.nochim, file=paste0(args[6],'dada2_asv_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  seq.out <- DNAStringSet(x=seqs, start=NA, end=NA, width=NA, use.names=TRUE)
  writeXStringSet(seq.out, filepath=paste(args[6],'representative_seq_set.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c( "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  write.table(track, file=paste0(args[6],'stats_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  seqtab2 <- seqtab.nochim
  #colnames(seqtab.nochim) == rownames(taxa) ! before changing -->  colnames(seqtab.nochim)<-new_names
 # taxa <- assignTaxonomy(seqtab.nochim, "/export/data01/databases/silva/r132/dada2/silva_nr_v132_train_set.fa.gz", multithread=15)
#  taxa2 <- addSpecies(taxa, "/export/data01/databases/silva/r132/dada2/silva_species_assignment_v132.fa.gz")
#  rownames(taxa2) <- new_names
}else{
  new_names <- c(paste("asv",1:length(colnames(seqtab2)),sep=""))
  seqs <- getSequences(seqtab2)
  #names(seqs) <- seqs
  names(seqs) <- new_names
  original_names<-colnames(seqtab2)
  colnames(seqtab2)<-new_names
  write.table(seqtab2, file=paste(args[6],'dada2_asv_table.txt',sep=""), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  seq.out <- DNAStringSet(x=seqs, start=NA, end=NA, width=NA, use.names=TRUE)
  writeXStringSet(seq.out, filepath=paste(args[6],'representative_seq_set.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2))
  colnames(track) <- c( "denoisedF", "denoisedR", "merged", "variants")
  rownames(track) <- sample.names
  write.table(track, file=paste0(args[6],'stats_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  
#  taxa <- assignTaxonomy(seqtab2, "/export/data01/databases/silva/r132/dada2/silva_nr_v132_train_set.fa.gz", multithread=15)
#  taxa2 <- addSpecies(taxa, "/export/data01/databases/silva/r132/dada2/silva_species_assignment_v132.fa.gz")
#  rownames(taxa2) <- new_names
  
  }
#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)

#TAXONOMY ASSIGNATION
extra_params_taxo <-args[14]
#if (!startsWith( trimws(extra_params_taxo), ',') && nchar(trimws(extra_params_taxo))>1){
#  extra_params <- paste(",",extra_params_taxo)
#}
#taxa <- eval(parse(text=paste("assignTaxonomy(seqtab2,args[11], multithread=cpus ",extra_params_taxo,")"))) 
colnames(seqtab2)<-original_names
taxa <- assignTaxonomy(seqtab2,args[11], multithread=cpus)
if (args[13] == "T" || args[13] == "TRUE" ){
  taxa2 <- addSpecies(taxa, args[12])
  rownames(taxa2) <- new_names
  write.table(taxa2, file=paste0(args[6],'taxonomy_dada2/representative_seq_set_tax_assignments.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  
}else{
  rownames(taxa) <- new_names
  write.table(taxa2, file=paste0(args[6],'taxonomy_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  
}

#to_assign <- readDNAStringSet("cascabel_project/runs/dada2/asv/representative_seq_set.fasta", format="fasta")
#newtax <- assignTaxonomy(to_assign,args[11], multithread=cpus)

#taxa22 <- addSpecies(newtax, args[12])
#rownames(taxa22) <- new_names
#write.table(taxa22, file=paste0(args[6],'taxonomy_dada2/representative_seq_set_tax_assignments_f.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)





