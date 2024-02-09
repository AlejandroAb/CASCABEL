                                                                    
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
#args[2] = cpus
#args[3] = plotErr
        #args[5] = extra parameters for dada2 function
#args[4] = path for outputfiles

#args[5]... = Number of errors
 #args[6]... = BinnedQ
#args[7]... = summary files from libraries

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
for(i in 7:(length(args))) {
  #paths <- c(paths,gsub('demultiplexed/','demultiplexed/filtered/',gsub("/summary.txt", '',args[i])))
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

#print(filtFs)

#if (args[2] == "pseudo"){
#  pool = "pseudo"
#}else{
#  pool <- eval(parse(text=args[2]))
#}
cpus <<- strtoi(args[2],10)
#nb <<- strtoi(args[5],10)
nb <- as.double(args[5])

print(nb)
print(args[5])
#extra_params <- args[5]

if (args[6] == "T" || args[6] == "TRUE" ){
   errF <- learnErrors(filtFs, multithread=cpus, nbases=nb, errorEstimationFunction = loessErrfun_mod1)
   errR <- learnErrors(filtRs, multithread=cpus, nbases=nb, errorEstimationFunction = loessErrfun_mod1)

}else{

  # learn error rates def.
  errF <- learnErrors(filtFs, multithread=cpus, nbases=nb)
  #errF <- eval(parse(text=paste("learnErrors(filtFs, multithread=cpus,",extra_params,")")))
  errR <- learnErrors(filtRs, multithread=cpus,nbases=nb)
}

#errR <- eval(parse(text=paste("learnErrors(filtRs, multithread=cpus)")))
#plotErrors(errF, nominalQ=TRUE)
if (args[3] == "T" || args[3] == "TRUE" ){
  library(ggplot2)
  plots_fw <- plotErrors(errF, nominalQ=TRUE)
  ggsave(paste(args[4],"fr_err.pdf",sep=""), plots_fw)
  plots_rv <- plotErrors(errR, nominalQ=TRUE)
  ggsave(paste(args[4],"rv_err.pdf", sep=""), plots_rv)
  #here we have to solve two things where to store the plots
  # and if we are going to generate one pdf per file, one per library
  # or one per strand...
}

errors_data <- paste(args[4],"errors.RData", sep="")
save(errF, errR, file = errors_data)
#exit 
quit("no",0)



