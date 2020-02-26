# Feb, 2020
# error correct barcodes with a levenstein distance up to args[4]
# write corrected fastq of the bc to file for qiime
# Script by J. Engelmann & A. Abdala

library(Biostrings)
library(ShortRead)
args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = file with barcode list,same as used for rule validate_mapping ($validate_mapping_file.py )
#args[3] = sbarcode.fq file -- result from extrac_barcodes.py
#args[4] = number of allowed mismatches
setwd(args[1])


# look up the expected barcodes in the sample sheet (they are the same for both libs)
ssheet   <- read.delim(paste('./',args[2], sep=''), stringsAsFactors=F)
exp.bars <- ssheet$BarcodeSequence
# rev complement of the barcodes, save as character string
#exp.barsRC <- as.character(reverseComplement(DNAStringSet(exp.bars)))
all.bars <- c(exp.bars)

mismatch = strtoi(args[4], 10)
### error correction ###
# read in barcodes with fastq streamer
# first compute edit (levenshtein) distances of observed barcodes to expected ones
# x is file name. this is somewhat slow ... but ok with nreads=5000
correctBC <- function(x, exp.barcodes, nreads){
	f <- FastqStreamer(x, nreads)
    	while(length(fq <- yield(f))) {
           dists <- adist(sread(fq), all.bars)
	   newBC <- as.vector(sread(fq))
	   for(i in 1:nrow(dists)){
	      if( min(dists[i,]) <= mismatch  &&  length(which(dists == min(dists)))<2){  # if min dist small enough and no ties, correct barcode (no change on perfect barcodes)
		  		newBC[i] <- exp.barcodes[which.min(dists[i,])]
		  # I cannot use sread to set values, need to construct a fqout, filling in sread, qual, whatever ...

		  }
	   }
	   newFQ      <- ShortReadQ(sread = DNAStringSet(newBC), quality = quality(fq), id = id(fq))
	   # mode='a': append, compress: should file be gzipped
	   writeFastq(newFQ, paste(x, "_corrected", sep=""), mode="a", compress=FALSE)
	 }
	 close(f)
}

correctBC(args[3], all.bars, 5000)

