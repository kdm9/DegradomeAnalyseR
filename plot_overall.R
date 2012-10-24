library(doMC)
# registerDoMC(cores=2)
library(foreach)
library(iterators)
source('~/Programming/Bioinformatics/DegradomeASC/DegradomeAnalyseR/peak_picker.R')

#setwd("/home/kevin/Downloads/ws/R")
args <- commandArgs(trailingOnly=T)
#args[1] = "deg_1_summarised.csv"
#args[2] = "deg_1_targets.csv"
#args[3] = "deg_1"
deg <- read.csv(args[1])
targets <- read.csv(args[2])

agis <- isplit(deg, deg$rname)

overall.hist <- data.frame(
    rname=character(),
    pos=numeric(),
    count=numeric(),
    m_mapq=numeric()
    )

foreach(agi=agis) %dopar% {
    seq.len = targets[targets[,1] == agi$key[[1]], 2]
    agi$value$pos <- agi$value$pos / seq.len
    overall.hist <- merge(overall.hist, data.frame(agi$value)[,2:5],all=T)
    write.csv(overall.hist, file=paste(args[3],"_normalised.csv", sep=""))
}
pdf(paste(args[3],"_overall.pdf", sep=""))
plot(overall.hist$pos, overall.hist$count)
dev.off()
write.csv(overall.hist, file=paste(args[3],"_normalised.csv", sep=""))

overall.peaks <- overall.hist
overall.peaks[with(overall.peaks, order(pos)), ]

peaks <- peak.pick(overall.peaks$count)
write.csv(peaks, file=paste(args[3],"_peaks.csv", sep=""))

overall.peaks$count <- peak.convert(overall.peaks$count)

sd(overall.hist$count)

write.csv(overall.peaks, file=paste(args[3],"_converted.csv", sep=""))

