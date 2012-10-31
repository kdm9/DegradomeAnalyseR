args <- commandArgs(trailingOnly=T)
# args[1] = "alx8_conserved.tab"
# args[2] = "alx8"

### READ IN DATA
# Read replicate cross-check filtered PAREsnip output
paresnip.filtered <- read.delim(args[1], header=T, sep='\t', quote="\"")

# Remove duplicate columns
SRA.FA.cols <- grep("^SRA|FA",names(paresnip.filtered))
paresnip.filtered <- paresnip.filtered[,c(1,2,3,SRA.FA.cols)]

# remove trailing ".1" from general column names 
names(paresnip.filtered)[1:4] <- gsub(
    "\\.1.*$",
    "",
    names(paresnip.filtered)[1:4]
    )    

# remove trailing whitespace from AGIs
paresnip.filtered$AGI <- gsub("\\s", "", paresnip.filtered$AGI)
paresnip.summarised$AGI <- gsub("\\s", "", paresnip.summarised$AGI)


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
