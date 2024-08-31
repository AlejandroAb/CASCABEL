
#####
# Dada2- learn errors and run dada2
#
# @author: J. Engelmann
# @author: A. Abdala
# @author: W. Reitsma

library(dada2)
library(Biostrings)


args <- commandArgs(trailingOnly = T)

#args[1] = path for seting WD to use current dir use $PWD
#args[2] = dada2 pool option
#args[3] = cpus
#args[4] = extra parameters for dada2 function
#args[5] = search chimera
#args[6] = path for outputfiles
#args[7] = shorts
#args[8] = longs
#args[9] = ofsets
#args[10] = minOverlap merge
#args[11] = maxMismatch merge
#args[12] = mergePairs
#args[13] = Interactive
#args[14] = NON Interactive behav.
#args[15] = chimera method
#args[16] = summary files from libraries



setwd(args[1])
paths <-NULL
for(i in 16:length(args)) {
  paths <- c(paths,gsub("summary.pcr.txt", 'filtered/',args[i]))
}
print(paths)
filtFs <-  sort(list.files(paths, pattern="_1.fastq.gz", full.names = TRUE))
filtRs <-  sort(list.files(paths, pattern="_2.fastq.gz", full.names = TRUE))
#Get sample names
sample.names <- gsub('_1.fastq.gz', '', basename(filtFs))

#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names


if (args[2] == "pseudo"){
  pool = "pseudo"
}else{
  pool <- eval(parse(text=args[2]))
}
cpus <<- strtoi(args[3],10)

# learn error rates
errors_data <- paste(args[6],"errors.RData", sep="")
load(errors_data)

extra_params <- args[4]
if (!startsWith( trimws(extra_params), ',') && nchar(trimws(extra_params))>1){
  extra_params <- paste(",",extra_params)
}

#merge parameters
minOv <- as.integer(args[10])
maxMism <- as.integer(args[11])
mergers <- vector("list")
dadaFs <- vector("list")
dadaRs <- vector("list")
for (sam in sample.names){
  cat("Processing", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose=TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose=TRUE)
  #dadaF <- eval(parse(text=paste("dada(derepF, err=errF, multithread=cpus, pool=FALSE, errorEstimationFunction = loessErrfun_mod1, ",extra_params,")")))   
  #dadaR <- eval(parse(text=paste("dada(derepR, err=errR, multithread=cpus, pool=FALSE, errorEstimationFunction = loessErrfun_mod1,",extra_params,")")))
  dadaF <- eval(parse(text=paste("dada(derepF, err=errF, multithread=cpus, pool=FALSE, ",extra_params,")")))
  dadaR <- eval(parse(text=paste("dada(derepR, err=errR, multithread=cpus, pool=FALSE, ", extra_params,")")))
  if (args[12] == "T" || args[12] == "TRUE" ){
    merger <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = minOv, maxMismatch = maxMism, returnRejects = TRUE, justConcatenate = TRUE)
  }else{
   merger <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = minOv, maxMismatch = maxMism)
  }
  mergers[[sam]] <- merger
  dadaFs[[sam]] <- dadaF
  dadaRs[[sam]] <- dadaR
}


# ASV table 
seqtab <- makeSequenceTable(mergers)

#rownames= samples #colnames= sequences

#this show sequence length distributions
seq_hist <- table(nchar(getSequences(seqtab)))
fname_seqh <- paste(args[6],"seq_hist.txt",sep="")
write.table(seq_hist, file = fname_seqh  , sep = "\t", quote=FALSE, col.names = FALSE)

shorts <- as.integer(args[7])
longs <- as.integer(args[8])
offset <- as.integer(args[9])
m <- mean(nchar(getSequences(seqtab)))
mx <- max(nchar(getSequences(seqtab)))
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
  cat("\n---This step can filter ASV based on user's expected fragmnt length---\n")
  cat("Please enter the option which fits better for your data:\n")
  cat(paste("1. Use values from configuration file: length > ", shorts," and length < ", longs,"\n"))
  cat(paste("2. Use values from median + /- the offset: length > ", (m-offset)," and length < ", (m+offset),"\n"))
  cat("3. Specify new values\n")
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
    
    shorts <<- readintegerConsole("Enter the shortest length allowed:")
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
    
    shorts <<- readinteger("Enter the shortest length allowed:")
    longs <<- readinteger("Enter the longest length allowed:")
  }
  else if(n == 4){
    print(seq_hist)
    return(createMenu())
  }else if(n == 5){
    quit("no",1)
  }else {
    return(createMenu())
  }
}

if (args[13] == "T" || args[13] == "TRUE" ){
  createMenuConsole()
}else if (args[14] == "AVG") {
  shorts <<- as.integer((m-offset))
  longs <<- as.integer((m+offset))
}else if (args[14] == "CFG") {
  #shorts <- strtoi(args[3], 10)
  #longs <- strtoi(args[4], 10)
  shorts <<- as.integer(args[7])
  longs <<- as.integer(args[8])
}else if (args[14] == "NONE") {
  #shorts <- strtoi(args[3], 10)
  #longs <- strtoi(args[4], 10)
  shorts <<- 0
  longs <<- mx + 1
}else{
 write("Invalid option for [rm_reads][non_interactive_behaviour] values at --configfile", stderr())
 quit("no",1)
}



#createMenuConsole()
write.table(c(shorts,longs),paste0(args[6],"shorts_longs.log"), quote=F,sep="\t", row.names = F, col.names = F)
#This is analogous to “cutting a band” in-silico to get amplicons of the targeted length
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% shorts:longs] 

if (args[5] == "T" || args[5] == "TRUE" ){
  #Pool by sample
  if (args[15] == "Consensus" || args[5] == "consensus" ){
    seqtab.nochim.tab <- isBimeraDenovoTable(seqtab2,  multithread=cpus, verbose=TRUE)
  }else{#Pool deNovo (this tends to generate more bimera)
    seqtab.nochim.tab <- isBimeraDenovo(seqtab2,  multithread=cpus, verbose=TRUE)
      
  }
# rm chimeras
  #seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=cpus, verbose=TRUE)  
  seqtab.nochim <- seqtab2[,names(seqtab.nochim.tab[seqtab.nochim.tab %in% c(F)])]
  seqtab.chim <- seqtab2[,names(seqtab.nochim.tab[seqtab.nochim.tab %in% c(T)])]
  #dim(seqtab.nochim)
  #sum(seqtab.nochim)/sum(seqtab)
  #matrix(1:9, nrow = 3, ncol = 3)
  seqs <- getSequences(seqtab.nochim)
  seqs.chim <- getSequences(seqtab.chim)
  #new_names <- c(paste("asv.",1:length(colnames(seqtab.nochim)),sep=""))
  new_names <- c(paste("asv.",1:length(seqs),sep=""))
  new_names.chim <- c(paste("bimera.",1:length(seqs.chim),sep=""))
  #estos son equivalentes colnames(seqtab.nochim) == seqs
  #names(seqs) <- seqs
  original_names<-seqs
  names(seqs) <- new_names
  original_names.chim<-seqs.chim
  #names(seqs.chim) <- new_names.chim
  #original_names<-colnames(seqtab.nochim)
  getN <- function(x) sum(getUniques(x))
  
  if(is.null(dim(seqtab.nochim))){
    tmp_seqtab.nochim <- matrix(as.list(seqtab.nochim),nrow=1)
    colnames(tmp_seqtab.nochim)<-new_names
    rownames(tmp_seqtab.nochim)<-sample.names
    write.table(tmp_seqtab.nochim, file=paste0(args[6],'dada2_asv_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
    
    #The normal output is a list with named elements, but working
    #with only one sample, this is no more a list and then we have
    #errors, following we re create the lists:
    tmp_dadaFs <- list(dadaFs)
    names(tmp_dadaFs) <- sample.names 
    tmp_dadaRs <- list(dadaRs)
    names(tmp_dadaRs) <- sample.names 
    tmp_mergers <- list(mergers)
    names(tmp_mergers) <- sample.names
    tmp_seqtab2 <- matrix(seqtab2,nrow=1)
    row.names(tmp_seqtab2) <-  sample.names
    tmp_seqtab.nochim2 <- matrix(seqtab.nochim,nrow=1)
    rownames(tmp_seqtab.nochim2)<-sample.names
    track <- cbind(sapply(tmp_dadaFs, getN), sapply(tmp_dadaRs, getN), sapply(tmp_mergers, getN), rowSums(tmp_seqtab2), rowSums(tmp_seqtab.nochim2))
    seqtab2 <- tmp_seqtab.nochim2
    
  }else{
    colnames(seqtab.nochim)<-new_names
    write.table(seqtab.nochim, file=paste0(args[6],'dada2_asv_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
    track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
    seqtab2 <- seqtab.nochim
    
  }
  if(dim(seqtab.chim)[2]==0){
    #write empty bimera table
    #colnames(seqtab.chim)<-new_names.chim
    #write.table(seqtab.chim, file=paste0(args[6],'dada2_asv_bimera_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
    file.create(paste(args[6],'dada2_asv_bimera_table.txt',sep=""))

  }else{
    print(dim(seqtab.chim))
    #write bimera, but we are not cheking for dimmension
    names(seqs.chim) <- new_names.chim

    colnames(seqtab.chim)<-new_names.chim
    write.table(seqtab.chim, file=paste0(args[6],'dada2_asv_bimera_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
    
  }
  #write bimera, but we are not cheking for dimmension
  #colnames(seqtab.chim)<-new_names.chim
  #write.table(seqtab.chim, file=paste0(args[6],'dada2_asv_bimera_table.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  # also write to file seqs
  seq.out <- DNAStringSet(x=seqs, start=NA, end=NA, width=NA, use.names=TRUE)
  seq.out.chim <- DNAStringSet(x=seqs.chim, start=NA, end=NA, width=NA, use.names=TRUE)
  
  writeXStringSet(seq.out, filepath=paste(args[6],'representative_seq_set.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  writeXStringSet(seq.out.chim, filepath=paste(args[6],'representative_bimeras.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  
  colnames(track) <- c( "denoisedF", "denoisedR", "merged","length", "nonchim")
  rownames(track) <- sample.names
  write.table(track, file=paste0(args[6],'stats_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
}else{
  new_names <- c(paste("asv",1:length(colnames(seqtab2)),sep=""))
  seqs <- getSequences(seqtab2)
  names(seqs) <- new_names
  original_names<-colnames(seqtab2)
  colnames(seqtab2)<-new_names
  write.table(seqtab2, file=paste(args[6],'dada2_asv_table.txt',sep=""), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  file.create(paste(args[6],'dada2_asv_bimera_table.txt',sep=""))
  seq.out <- DNAStringSet(x=seqs, start=NA, end=NA, width=NA, use.names=TRUE)
  writeXStringSet(seq.out, filepath=paste(args[6],'representative_seq_set.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  #create empty file
  file.create(paste(args[6],'representative_bimeras.fasta',sep=""))

  getN <- function(x) sum(getUniques(x))
  
  ####
  if(is.null(dim(seqtab2))){
    tmp_dadaFs <- list(dadaFs)
    names(tmp_dadaFs) <- sample.names 
    tmp_dadaRs <- list(dadaRs)
    names(tmp_dadaRs) <- sample.names 
    tmp_mergers <- list(mergers)
    names(tmp_mergers) <- sample.names
    tmp_seqtab2 <- matrix(seqtab2,nrow=1)
    row.names(tmp_seqtab2) <-  sample.names
  track <- cbind(sapply(tmp_dadaFs, getN), sapply(tmp_dadaRs, getN), sapply(tmp_mergers, getN), rowSums(tmp_seqtab2))
  }else{
    
    track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab2))
    
  }
  ###
  
  colnames(track) <- c( "denoisedF", "denoisedR", "merged", "variants")
  rownames(track) <- sample.names
  write.table(track, file=paste0(args[6],'stats_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  

}




