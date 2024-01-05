
#####
# Dada2- learn errors and run dada2
#
# @author: J. Engelmann
# @author: A. Abdala
# @author: W. Reitsma

library(dada2)
library(Biostrings)
library(magrittr)
library(dplyr)


args <- commandArgs(trailingOnly = T)

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
#args[14]... = extra_params taxonomy
#args[15]... = minOverlap merge
#args[16]... = maxMismatch merge
#args[17]... = mergePairs
#args[18]... = add_sps extra_params
#args[19]... = Interactive
#args[20]... = NON Interactive behav.
#args[21]... = summary files from libraries

# error function model obtained from Hannah Holland-Moritz on dada2 github:
# https://github.com/benjjneb/dada2/issues/1307#issuecomment-957680971
loessErrfun_mod1 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}


setwd(args[1])
paths <-NULL
for(i in 21:(length(args)-1)) {
  #paths <- c(paths,gsub('demultiplexed/','demultiplexed/filtered/',gsub("/summary.txt", '',args[i])))
   paths <- c(paths,gsub("summary.pcr.txt", 'filtered/',args[i]))
}
print(paths)
filtFs <-  sort(list.files(paths, pattern="_1.fastq.gz", full.names = TRUE))
filtRs <-  sort(list.files(paths, pattern="_2.fastq.gz", full.names = TRUE))
#Get sample names
sample.names <- gsub('_1.fastq.gz', '', basename(filtFs))
print(sample.names)
#assign names to files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

print(filtFs[1])

if (args[2] == "pseudo"){
  pool = "pseudo"
}else{
  pool <- eval(parse(text=args[2]))
}
cpus <<- strtoi(args[3],10)
extra_params <- args[5]

# learn error rates
errF <- learnErrors(filtFs, multithread=cpus, nbases=1e8, errorEstimationFunction = loessErrfun_mod1)
errR <- learnErrors(filtRs, multithread=cpus, nbases=1e8, errorEstimationFunction = loessErrfun_mod1)

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
if (!startsWith( trimws(extra_params), ',') && nchar(trimws(extra_params))>1){
  extra_params <- paste(",",extra_params)
}

#merge parameters
minOv <- as.integer(args[15])
maxMism <- as.integer(args[16])
mergers <- vector("list")
dadaFs <- vector("list")
dadaRs <- vector("list")
for (sam in sample.names){
  cat("Processing", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]], verbose=TRUE)
  derepR <- derepFastq(filtRs[[sam]], verbose=TRUE)
  dadaF <- eval(parse(text=paste("dada(derepF, err=errF, multithread=cpus, pool=FALSE,  errorEstimationFunction = loessErrfun_mod1, ",extra_params,")")))   
  dadaR <- eval(parse(text=paste("dada(derepR, err=errR, multithread=cpus, pool=FALSE, errorEstimationFunction = loessErrfun_mod1,",extra_params,")")))
  if (args[17] == "T" || args[17] == "TRUE" ){
    merger <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = minOv, maxMismatch = maxMism, returnRejects = TRUE, justConcatenate = TRUE)
  }else{
   merger <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap = minOv, maxMismatch = maxMism)
  }
  mergers[[sam]] <- merger
  dadaFs[[sam]] <- dadaF
  dadaRs[[sam]] <- dadaR
}


print(names(dadaFs))

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

if (args[19] == "T" || args[19] == "TRUE" ){
  createMenuConsole()
}else if (args[20] == "AVG") {
  shorts <<- as.integer((m-offset))
  longs <<- as.integer((m+offset))
}else if (args[20] == "CFG") {
  #shorts <- strtoi(args[3], 10)
  #longs <- strtoi(args[4], 10)
  shorts <<- as.integer(args[7])
  longs <<- as.integer(args[8])
}else if (args[20] == "NONE") {
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

if (args[10] == "T" || args[10] == "TRUE" ){
   
# rm chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=cpus, verbose=TRUE)  
  #dim(seqtab.nochim)
  #sum(seqtab.nochim)/sum(seqtab)
  #matrix(1:9, nrow = 3, ncol = 3)
  seqs <- getSequences(seqtab.nochim)
  #new_names <- c(paste("asv.",1:length(colnames(seqtab.nochim)),sep=""))
  new_names <- c(paste("asv.",1:length(seqs),sep=""))
  #estos son equivalentes colnames(seqtab.nochim) == seqs
  #names(seqs) <- seqs
  original_names<-seqs
  names(seqs) <- new_names
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
  # also write to file
  seq.out <- DNAStringSet(x=seqs, start=NA, end=NA, width=NA, use.names=TRUE)
  writeXStringSet(seq.out, filepath=paste(args[6],'representative_seq_set.fasta',sep=""), append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
  colnames(track) <- c( "denoisedF", "denoisedR", "merged","length", "nonchim")
  rownames(track) <- sample.names
  write.table(track, file=paste0(args[6],'stats_dada2.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
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
  
#  taxa <- assignTaxonomy(seqtab2, "/export/data01/databases/silva/r132/dada2/silva_nr_v132_train_set.fa.gz", multithread=15)
#  taxa2 <- addSpecies(taxa, "/export/data01/databases/silva/r132/dada2/silva_species_assignment_v132.fa.gz")
#  rownames(taxa2) <- new_names
  
}
#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)

#TAXONOMY ASSIGNATION
extra_params_taxo <-args[14]
if (!startsWith( trimws(extra_params_taxo), ',') && nchar(trimws(extra_params_taxo))>1){
  extra_params_taxo <- paste(",",extra_params_taxo)
}
#dadaFs <- eval(parse(text=paste("dada(filtFs, err=errF, multithread=cpus, pool=pool, ",extra_params,")")))
colnames(seqtab2)<-original_names
taxa <- eval(parse(text=paste("assignTaxonomy(seqtab2,args[11], multithread=cpus ",extra_params_taxo,")"))) 
#taxa <- assignTaxonomy(seqtab2,args[11], multithread=cpus)
if (args[13] == "T" || args[13] == "TRUE" ){
 # taxa2 <- addSpecies(taxa, args[12])
  extra_params_add_sp <-args[18]
  if (!startsWith( trimws(extra_params_add_sp), ',') && nchar(trimws(extra_params_add_sp))>1){
    extra_params_add_sp <- paste(",",extra_params_add_sp)
  }
  taxa2 <- eval(parse(text=paste("addSpecies(taxa, args[12]", extra_params_add_sp,")")))
  rownames(taxa2) <- new_names
  write.table(taxa2, file=paste0(args[6],'taxonomy_dada2/representative_seq_set_tax_assignments.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  
}else{
  rownames(taxa) <- new_names
  write.table(taxa, file=paste0(args[6],'taxonomy_dada2/representative_seq_set_tax_assignments.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  
}

#to_assign <- readDNAStringSet("cascabel_project/runs/dada2/asv/representative_seq_set.fasta", format="fasta")
#newtax <- assignTaxonomy(to_assign,args[11], multithread=cpus)

#taxa22 <- addSpecies(newtax, args[12])
#rownames(taxa22) <- new_names
#write.table(taxa22, file=paste0(args[6],'taxonomy_dada2/representative_seq_set_tax_assignments_f.txt'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)





