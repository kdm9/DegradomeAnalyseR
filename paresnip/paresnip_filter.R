# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the
# GPLv3 License

library(plyr)
args <- commandArgs(trailingOnly=T)
#args should be the name of each file, followed by the genotype of the files
# args[1] = "paresnip_1.tab"
# args[2] = "paresnip_2.tab"
# args[3] = "paresnip_2.tab"
# args[4] = "col0"

cross.validate <- function (params){
  counter = 0
  summed.data = list()
  replicated.data = list()

  # This first loop loads the replicates into a single data matrix
  for (param in params){
    counter = counter + 1
    data <- read.delim(param[[1]],header=T, quote = '"', sep='\t')

    # Exract AGIs from fasta header lines which include descriptions
    agi.list <- unlist(strsplit(as.character(data$Gene), "|", fixed=T))
    data$AGI <- matrix(agi.list[grep("^AT.G", agi.list)])

    # Collapse data
    data.summarised <- ddply(
      data[data$Category<3,],
      .(AGI, Cleavage.Position, Duplex),
      summarize,
      FA=max(Fragment.Abundance),
      SRA=max(Short.Read.Abundance)
    )

    # Create key, this should be the lowest common denominator of the
    # collapsed data, and therefore should be unique
    data.summarised$Key <- paste(
      data.summarised$AGI,
      data.summarised$Cleavage.Position,
      data.summarised$Duplex
      )

    # Sort by keys, so that the match() below works
    data.summarised <- data.summarised[with(data.summarised, order(Key)), ]

    # Fix names
    names(data.summarised) <- paste(names(data.summarised), counter)

    # Append to data matrix
    summed.data <- c(summed.data, data.summarised)
    print(length(summed.data))

  }

  print(names(summed.data))

  # This second loop does the match operation, keeping only data which is
  # replicated across the samples
  for (iii in c(1:(counter))){
    cols = 6
    keycol = 6

    # Add the first bit of the big data matrix to the replicated data matrix
    # and skip to the next iteration, if this is the first iteration of the
    # loop
    if (iii == 1){
      replicated.data <- data.frame(summed.data[1:cols])
      next
    }

    # Take a chunk from the big data matrix
    data <- data.frame(
        summed.data[seq((iii*cols)-(cols-1),(iii*cols))]
        )

    # Remove rows from the chunk of data which are not present in the
    # replicated data
    a = as.character(replicated.data[,length(replicated.data)])
    b = as.character(data$Key)
    matches <- match(a, b)
    matches <- matches[!is.na(matches)]
    data <- data[matches,]

    # Remove rows from replicated data not present in the chunk
    a = as.character(replicated.data[,length(replicated.data)])
    b = as.character(data$Key)
    matches <- match(b, a)
    matches <- matches[!is.na(matches)]
    replicated.data <- replicated.data[matches,]

    # Append chunk to the replicated data matrix
    replicated.data <- data.frame(replicated.data, data)
  }
  return(replicated.data)
}

# Grab file list from commandline
file.list <- list()
for (file in args[seq(1,(length(args)-1))]){
  file.tuple <- list(file,args[length(args)])
  file.list <- c(file.list, list(file.tuple))
}

fn = paste(args[length(args)], "conserved.tab", sep=".")
conserved <- data.frame(cross.validate(file.list))
write.table(
  conserved,
  file=fn,
  sep="\t",
  col.names=T,
  row.names=F, # Don't write row names
  qmethod = "double" # quote method
)
