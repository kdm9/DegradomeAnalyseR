library(doMC)
# registerDoMC(cores=2)
library(foreach)
library(iterators)

args <- commandArgs(trailingOnly=T)
# args[1] = "alx8_summarised.csv"
# args[2] = "header.csv"
# args[3] = "alx8_rnaseq.filtered.csv"
# args[4] = "alx8"

### READ DATA
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


### CALCULATE BIAS TABLE
# pre-fill summarised table based on all agis
RNAseq.cleavage.summary <- data.frame(
    agi=agi.lengths$agi,
    upstream.count=numeric(length(agi.lengths$agi)),
    downstream.count=numeric(length(agi.lengths$agi)),
    upstream.length=numeric(length(agi.lengths$agi)),
    downstream.length=numeric(length(agi.lengths$agi)),
    )

RNAseq.cleavage.bias <- data.frame(
  agi=agi.lengths$agi,
  bias=numeric(length(agi.lengths$agi)),
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
    
    # Fill RNAseq cleavage summary table
    RNAseq.cleavage.summary[RNAseq.cleavage.summary[,1]==agi,] <- c(
        agi,
        upstream.count,
        downstream.count,
        upstream.length,
        downstream.length
        )
    
    # calculate bias
    upstream.RPK <- upstream.count / (upstream.length / 1000)
    downstream.RPK <- downstream.count / (downstream.length / 1000)
    
    bias <- downstream.RPK / upstream.RPK
    
    # Fill RNAseq cleavage bias table
    RNAseq.cleavage.bias[RNAseq.cleavage.bias[,1]==agi,] <- c(
      agi,
      bias
    )
    
    
}

write.csv(RNAseq.cleavage.summary,file=paste(args[4], "cleavage.bias.csv", sep="."))
write.csv(RNAseq.cleavage.bias,file=paste(args[4], "cleavage.bias.csv", sep="."))

## DO TESTS
