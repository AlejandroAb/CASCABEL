#####
# Dada2 - assign taxonomy
#
# @author: J. Engelmann
# @author: A. Abdala

library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = T)

#args[1]= "/export/lv11/projects/Marco-Bolo/Data/Operational/CASCABEL/"   #Working directory
#args[2]= "/export/lv11/projects/Marco-Bolo/Data/Operational/CASCABEL/Marco-bolo/runs/16S_max_overlap_l200/asv/representative_bimeras.fasta"
#args[3] ="minBoot=60"
#args[4] = "allowMultiple=TRUE"
#args[5] = 4249
#args[6] = "/export/lv11/projects/Marco-Bolo/Data/RawData/16S/Reference_library_16S_AllFish_Dez2023_LineageInfo_sintax_format_6levels_dad2.fasta"
#args[7] = "/export/lv11/projects/Marco-Bolo/Data/RawData/16S/Reference_library_16S_AllFish_Dez2023_LineageInfo_sintax_format_sps_dad2.fasta"
#args[8] = "T"
#args[9] = "Marco-bolo/runs/16S_max_overlap_l200/asv/taxonomy_dada2/bimeras_tax_assignments.txt"
#args[10] = "Marco-bolo/runs/16S_max_overlap_l200/asv/taxonomy_dada2/bimeras_tax_assignments.bootstrap.txt"
#args[11]= 20
#args[12]= "bimera."
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
#new_names <- c(paste("asv.",1:length(s),sep=""))
new_names <- c(paste(args[12],1:length(s),sep=""))

conf.seed <<- strtoi(args[5],10)

set.seed(conf.seed)
if(length(seqFromFile)>0){
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
    rownames(taxa$tax) <- new_names
    write.table(taxa$tax, file=args[9], sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  }
}else{
  file.create(args[9])
  file.create(args[10])
}

