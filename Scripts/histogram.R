#FUnction to compute the Mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#Read arguments from command line.
args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD
#args[2] = file with sequence distribution
#args[3] = file with sequence lengths
#args[4] = path to store the files
#args[5] = Name of the file seqs_dist_hist.{sample}.png

setwd(args[1])

file <- args[2]
histTable <- read.table(file)

file2 <- args[3]
lengthTable <- read.table(file2)
summ <-summary(lengthTable$V1)
#summ <- summary(histTable$V2)
moda <- Mode(lengthTable$V1)
summ[7] <- moda
write(summ, paste(args[4],"seqs_statistics.txt", sep=''), ncolumns = 7)

#png(paste(args[4],"seqs_dist_hist.png", sep=''))
png(args[5])
plot(histTable$V2, histTable$V1, type="l", main="Sequence length distribution",xlab="Length", ylab="Number of sequences")
#dev.off()
