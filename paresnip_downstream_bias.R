library(doMC)
# registerDoMC(cores=2)
library(foreach)
library(iterators)

args <- commandArgs(trailingOnly=T)
# args[1] = "col0.summarised.csv"
# args[2] = "tx_tair10rep.targets.csv"
# args[3] = "rnaseq_1.filtered.csv"
# args[4] = "col0"

### READ DATA
agi.peaks.raw <- read.csv(args[1])
agi.peaks <- agi.peaks.raw[,-(1)]
row.names(agi.peaks) <- agi.peaks.raw$X
#[,1]: agi
#[,2]: peak pos or NA
#[,3]: paresnip category
# remaining columns are quantification data

agi.lengths <- read.csv(args[2])
names(agi.lengths) <- c("agi", "length")
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
    agi=as.character(agi.lengths$agi),
    upstream.count=numeric(length(agi.lengths$agi)),
    downstream.count=numeric(length(agi.lengths$agi)),
    upstream.length=numeric(length(agi.lengths$agi)),
    downstream.length=numeric(length(agi.lengths$agi))
    )

RNAseq.cleavage.bias <- data.frame(
  agi=agi.lengths$agi,
  bias=numeric(length(agi.lengths$agi))
)

# for every agi in the rnaseq bam, find the number of reads upstream and 
# downstream of the cleavage site, then find RPK of these
RNAseq.iter <- isplit(RNAseq.filtered, RNAseq.filtered$rname)

foreach (RNAseq.agi=RNAseq.iter) %dopar%{
    
    agi <- RNAseq.agi$key[[1]]
    agi.length <- agi.lengths$length[agi.lengths$agi==agi]
    
    # get a peak position, if there's no peak, get random position
    peak.position <- mean(agi.peaks[agi.peaks[,1]==agi, 2])
    if (is.na(peak.position)){
        peak.position <- runif(1) * agi.length
    }
    if (is.na(peak.position)){
        next
    }
    
    upstream.length <- peak.position
    downstream.length <- agi.length - peak.position
    
    positions <- as.numeric(RNAseq.agi$value$pos)
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

# Remove 0's and INF's, as these mean data is missing
RNAseq.cleavage.bias$bias[RNAseq.cleavage.bias$bias==0] <- NA
RNAseq.cleavage.bias$bias[RNAseq.cleavage.bias$bias==Inf] <- NA

#write data
write.csv(RNAseq.cleavage.summary,file=paste(args[4], "cleavage.bias.csv", sep="."))
write.csv(RNAseq.cleavage.bias,file=paste(args[4], "cleavage.bias.csv", sep="."))

## DO TESTS
RNAseq.cleavage.bias <- RNAseq.cleavage.bias[!is.na(RNAseq.cleavage.bias[,2]),]
table(!is.na(RNAseq.cleavage.bias$bias))
peak.agis <- unique(as.character(agi.peaks$AGI))
peak.bias.agi.match <- match(as.character(RNAseq.cleavage.bias$agi),peak.agis)
RNAseq.cleavage.bias$bias[peak.bias.agi.match] = RNAseq.cleavage.bias$bias[peak.bias.agi.match] -1.9
peak.bias.agi.match <- peak.bias.agi.match[!is.na(peak.bias.agi.match)]
peak.biases <- as.numeric(RNAseq.cleavage.bias$bias[peak.bias.agi.match])
nopeak.biases <- sample(as.numeric(RNAseq.cleavage.bias$bias[-(peak.bias.agi.match)]),length(peak.biases))

hist(peak.biases, breaks=100)
hist(nopeak.biases, breaks=100)

t.test(peak.biases, nopeak.biases)

wilcox.test(peak.biases, nopeak.biases)