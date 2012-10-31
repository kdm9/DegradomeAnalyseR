library(doMC)
# registerDoMC(cores=2)
library(foreach)
library(iterators)

args <- commandArgs(trailingOnly=T)

## READ DATA

agi.peaks <- read.csv(args[1])
#[,1]: agi
#[,2]: peak pos or NA
#[,3]: paresnip category
# remaining columns are quantification data

agi.lengths <- read.csv(args[2])
#[,1]: agi
#[,2]: agi length

RNAseq.filtered <- read.csv(args[3])
#[,1]: rname, which is an agi
#[,2]: pos, mapping position of read
#[,3]: mapq, mapping qualtity score
#[,4]: cigar, mapping indel string


## CALCULATE BIAS TABLE
# pre-fill summarised table based on all agis
RNAseq.cleavage.bias <- data.frame(
    agi=agi.lengths$agi,
    upstream.positions=numeric(length(agi.lengths$agi)),,
    downstream.positions=numeric(length(agi.lengths$agi)),,
    upstream.length=numeric(length(agi.lengths$agi)),,
    downstream.length=numeric(length(agi.lengths$agi)),
    )

# for every agi in the rnaseq bam, find the number of reads upstream and 
# downstream of the cleavage site, then find RPK of these
RNAseq.iter <- isplit(deg, deg$rname)
foreach (RNAseq.agi=RNAseq.iter) %dopar%{
    
    agi <- RNAseq.agi$key[[1]]
    agi.length <- agi.lengths[agi.lengths[,1]==agi, 2]
    
    # get a peak position, if there's no peak, get random position
    peak.position <- agi.peaks[agi.peaks[,1]==agi, 2]
    if (is.na(peak.position)){
        peak.position <- runif(1) * agi.length
    }
    if (is.na(peak.position)){
        next
    }
    
    upstream.length <- peak.position
    downstream.length <- agi.length - peak.position
    
    positions <- as.numeric(RNAseqagi$value$pos)
    upstream.count <- length(positions[positions<peak.position])
    downstream.count <- length(positions[positions>=peak.position])
    
    # Fill RNAseq cleavage bias table
    RNAseq.cleavage.bias[RNAseq.cleavage.bias[,1]==agi,] <- c(
        agi,
        upstream.positions,
        downstream.positions,
        upstream.length,
        downstream.length
        )
}

write.csv(RNAseq.cleavage.bias,paste(args[4], "cleavage_bias", sep="_"))

## DO TESTS

