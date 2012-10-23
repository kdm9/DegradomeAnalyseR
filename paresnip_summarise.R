args <- commandArgs(trailingOnly=T)

paresnip.filtered <- read.delim(args[1], header=T, sep='\t', quote="\"")
SRA.FA.cols <- grep("^SRA|FA",names(paresnip.filtered))
paresnip.filtered <- paresnip.filtered[,c(1,2,3,SRA.FA.cols)]

paresnip.summarised <- paresnip.filtered
SRA.cols <- grep("^SRA",names(paresnip.summarised))
FA.cols <- grep("^FA",names(paresnip.summarised))
paresnip.summarised$SRA.mean <- apply(paresnip.summarised[,SRA.cols], 1, mean)
paresnip.summarised$FA.mean <- apply(paresnip.summarised[,FA.cols], 1, mean)

paresnip.summarised$SRA.sd <- apply(paresnip.summarised[,SRA.cols], 1, sd)
paresnip.summarised$FA.sd <- apply(paresnip.summarised[,FA.cols], 1, sd)

write.csv(paresnip.summarised, paste(args[2], "summarised.csv", sep="_"))