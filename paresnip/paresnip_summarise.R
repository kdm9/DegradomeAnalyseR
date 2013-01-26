# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the
# GPLv3 License

args <- commandArgs(trailingOnly=T)
# args[1] = "alx8_conserved.tab"
# args[2] = "alx8"

### READ IN DATA
# Read replicate cross-check filtered PAREsnip output
paresnip.filtered <- read.delim(args[1], header=T, sep='\t', quote="\"")

# Remove duplicate columns
SRA.FA.cols <- grep("^SRA|FA",names(paresnip.filtered))
# 1=AGI, 2=cleavage.position, 3=duplex
general.cols <- c(1,2,3)
paresnip.filtered <- paresnip.filtered[,c(general.cols,SRA.FA.cols)]

# remove trailing ".1" from general column names
names(paresnip.filtered)[general.cols] <- gsub(
    "\\.1.*$",
    "",
    names(paresnip.filtered)[general.cols]
    )

# remove trailing whitespace from AGIs
paresnip.filtered$AGI <- gsub("\\s", "", paresnip.filtered$AGI)


### SUMMARISE DATA
paresnip.summarised <- paresnip.filtered

# Get SRA and FA column numbers
SRA.cols <- grep("^SRA",names(paresnip.summarised))
FA.cols <- grep("^FA",names(paresnip.summarised))

# Calculate Means and SDs of SRA and FA replicates.
paresnip.summarised$SRA.mean <- apply(paresnip.summarised[,SRA.cols], 1, mean)
paresnip.summarised$FA.mean <- apply(paresnip.summarised[,FA.cols], 1, mean)
paresnip.summarised$SRA.sd <- apply(paresnip.summarised[,SRA.cols], 1, sd)
paresnip.summarised$FA.sd <- apply(paresnip.summarised[,FA.cols], 1, sd)


write.csv(paresnip.summarised, paste(args[2], "summarised.csv", sep="."))
