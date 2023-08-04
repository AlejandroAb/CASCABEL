#Read arguments from command line.
args <- commandArgs(trailingOnly = T)
#args[1] = path for seting WD to use current dir use $PWD 
#args[2] = file with matrix 
#args[3] = heatmap png out 
#args[4] = heatmap png out with golay names (if args[2]"_golay" exists)
setwd(args[1]) 
mtx_file <- args[2] 
bcs <- read.delim(mtx_file, header=TRUE, row.names=1) 
bcs_matrix=as.matrix(bcs)
png(args[3]) 
heatmap(bcs_matrix, Colv = NA, Rowv = NA,scale = "column", main="Number of reads per barcode pair",xlab="Forward barcode", ylab="Reverse barcode") 
graphics.off() 

golay_mtx <- paste(mtx_file,'_golay',sep='') 
if (file.exists(golay_mtx)){
  bcs_golay <- read.delim(golay_mtx, header=TRUE, row.names=1)
  golay_mtx<-as.matrix(bcs_golay)
  png(args[4])
  heatmap(golay_mtx,Colv = NA, Rowv = NA,scale = "column",main="Number of reads per barcode pair",xlab="Forward barcode", ylab="Reverse barcode")
  graphics.off()
  
}
