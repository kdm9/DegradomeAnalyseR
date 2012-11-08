# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the
# GPLv3 License

library(plyr)
args <- commandArgs(trailingOnly=T)

##Import filtered bam
bam <- read.csv(args[1])

## Summarise data to count table
bam.counts <- ddply(
    bam,
    .(rname, pos),
    summarise, 
    count=length(rname),
    m_mapq=mean(mapq)
    )

## Write filtered bam to CSV
csv_filename <- paste(args[2], "summarised.csv", sep=".")
write.csv(bam.counts, file=csv_filename)
