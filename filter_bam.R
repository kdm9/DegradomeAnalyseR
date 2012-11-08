# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the GPLv3 License

library(Rsamtools)
args <- commandArgs(trailingOnly=T)

### EXPORT BAM

## Import data from bam
#rname:	reference name
#pos:	ref. mapping position
#mapq:	mapping quality score
#cigar:	CIGAR indel string
bam <- scanBam(
    args[1], 
    param=ScanBamParam(what=c("rname","pos", "mapq", "cigar"))
    )
bam <- data.frame(bam)

## Filter data
#Remove NA rows (unmapped reads)
bam <- bam[!is.na(bam[,1]),]
bam <- bam[!is.na(bam[,2]),]
#Remove rows with mapq < 20
# bam <- bam[bam[,3]>=20,]

## Write filtered bam to CSV
csv_filename <- paste(args[2], "filtered.csv", sep=".")
write.csv(bam, file=csv_filename)


### EXPORT BAM TARGETS
## import header
header <- scanBamHeader(args[1])
#extract targets detail
targets <- data.frame(header[[1]]$targets)

targets_filename <- paste(args[2], "targets.csv", sep=".")
write.csv(targets, file=targets_filename)

