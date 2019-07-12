# Aug 9, 2017
# error correct barcodes with a levenstein distance up to 2
# write corrected fastq of the bc to file for qiime
# and demultiplexed read fastq for dada2.
# Script by J. Engelmann
#Parametrization by A. Abdala

library(Biostrings)
library(ShortRead)
args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = file with barcode list,same as used for rule validate_mapping ($validate_mapping_file.py )
#args[3] = sbarcode.fq file -- result from extrac_barcodes.py
#args[4] = number of allowed mismatches
#AA setwd("~/Projects/eDNA_NIOZ47_51/Data/Operational/BarcodesExtracted")
setwd(args[1])

# read in NIOZ47 and 51 barcodes
#fnames  <- c('NIOZ47/barcodes.fastq', 'NIOZ51/barcodes.fastq')
#libs    <- gsub('/barcodes.fastq', '', fnames)  # need lib names later

# look up the expected barcodes in the sample sheet (they are the same for both libs)
#AA ssheet   <- read.delim('~/Projects/eDNA_NIOZ47_51/sampleList_mergedBarcodes_NIOZ47.txt', stringsAsFactors=F)
ssheet   <- read.delim(paste('./',args[2], sep=''), stringsAsFactors=F)
exp.bars <- ssheet$BarcodeSequence
# rev complement of the barcodes, save as character string
exp.barsRC <- as.character(reverseComplement(DNAStringSet(exp.bars)))
all.bars <- c(exp.bars, exp.barsRC)

# (check that the minimum distance b/w any expected barcode is larger 3 ->
# don't need to do that bc same barcodes as NIOZ70-71)


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
	      if( min(dists[i,]) <= args[4]){  # if min dist small enough, correct barcode (no change on perfect barcodes)
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


# correct the barcodes of the two files. Then check how many reads got assigned.
#AAcorrectBC(fnames[1], all.bars, 5000)
#AAcorrectBC(fnames[2], all.bars, 5000)
#correctBC(paste(paste("data/barcodes/",args[3],sep=''),'/barcodes.fastq',sep=''), all.bars, 5000)
#correctBC(paste(paste("data/barcodes/",args[3],sep=''),'/barcodes.fastq',sep=''), all.bars, 5000)
correctBC(args[3], all.bars, 5000)
# mv the filenames to corrected_barcodes.fastq on cmd line.

# use fastq with corrected barcodes for split_libraries.py on merged files!


##### demultiplex ### use for dada2 ### not yet run
# writes one fq file for each barcode sequence. x is multiplexed file.
# does so in chunks of nreads.
demultiplex <- function(x, barcode, nreads) {
    f <- FastqStreamer(x, nreads)
    while(length(fq <- yield(f))) {
        for(i in barcode) {
           #  pattern <- paste("^", i, sep="")
           #  fqsub <- fq[grepl(pattern, sread(fq))]

	   fqsub <- fq[sread(fq)==barcode]
	   if(length(fqsub) > 0) writeFastq(fqsub, paste(x, i, sep="_"), mode="a")
        }
    }
    close(f)
}

# demultiplex(x=fastq[1], barcode=c("TT", "AA", "GG"), nreads=500)
