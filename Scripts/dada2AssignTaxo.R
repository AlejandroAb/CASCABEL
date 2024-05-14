#####
# Dada2 - assign taxonomy
#
# @author: J. Engelmann
# @author: A. Abdala

library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = T)

#args[1]= "/export/lv5/scratch/projects_AA/CASCABEL_6/CASCABEL"   #Working directory
#args[2]= Path to the fasta file for assignation
#args[3] = extra_params_taxo
#args[4] = extra_params_taxo_sp
#args[5] = seed.
#args[6] = db
#args[7] = db. sp
#args[8] = add sp?
#args[9] = out assignment
#args[10] = out bootstrap
#args[11]= cpus

#TAXONOMY ASSIGNATION

setwd(args[1])
cpus <<- strtoi(args[11],10)
extra_params_taxo <-args[3]
if (!startsWith( trimws(extra_params_taxo), ',') && nchar(trimws(extra_params_taxo))>1){
  extra_params_taxo <- paste(",",extra_params_taxo)
}
#dadaFs <- eval(parse(text=paste("dada(filtFs, err=errF, multithread=cpus, pool=pool, ",extra_params,")")))
#colnames(seqtab2)<-original_names
#taxa <- eval(parse(text=paste("assignTaxonomy(seqtab2,args[11], multithread=cpus ",extra_params_taxo,")")))
#taxa <- assignTaxonomy(seqtab2,args[11], multithread=cpus)
#seqFromFile <- readDNAStringSet(filepath=paste(args[6],'representative_seq_set.fasta',sep=""))
seqFromFile <- readDNAStringSet(filepath=args[2])

s <- getSequences(seqFromFile);
new_names <- c(paste("asv.",1:length(s),sep=""))

conf.seed <<- strtoi(args[5],10)

set.seed(conf.seed)

taxa <- eval(parse(text=paste("assignTaxonomy(s,args[6],outputBootstraps = TRUE, multithread=cpus ",extra_params_taxo,")")))
rownames(taxa$boot) <- new_names
write.table(taxa$boot, file=args[10], sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
if (args[8] == "T" || args[8] == "TRUE" ){
 # taxa2 <- addSpecies(taxa, args[12])
  extra_params_add_sp <-args[4]
  if (!startsWith( trimws(extra_params_add_sp), ',') && nchar(trimws(extra_params_add_sp))>1){
    extra_params_add_sp <- paste(",",extra_params_add_sp)
  }
  taxa2 <- eval(parse(text=paste("addSpecies(taxa$tax, args[7]", extra_params_add_sp,")")))
  rownames(taxa2) <- new_names
  write.table(taxa2, file=args[9], sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)

}else{
  rownames(taxa) <- new_names
  write.table(taxa$tax, file=args[9], sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
}

