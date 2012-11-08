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
  for (param in params){
    counter = counter + 1
    data <- read.delim(param[[1]],header=T, quote = '"', sep='\t')
    agi.list <- unlist(strsplit(as.character(data$Gene), "|", fixed=T))
    data$AGI <- matrix(agi.list[grep("^AT.G", agi.list)])

    data.summarised <- ddply(
      data[data$Category<3,],
      .(AGI, Cleavage.Position, Duplex),
      summarize,
      FA=max(Fragment.Abundance),
      SRA=max(Short.Read.Abundance),
    )

    data.summarised$Key <- paste(
      data.summarised$AGI,
      data.summarised$Cleavage.Position
      data.summarised$Duplex
      )

    data.summarised <- data.summarised[with(data.summarised, order(Key)), ]
    names(data.summarised) <- paste(names(data.summarised), counter)
    summed.data <- c(summed.data, data.summarised)
    print(length(summed.data))

  }
  print(names(summed.data))
  for (iii in c(1:(counter))){
    cols = 6
    keycol = 7
    if (iii == 1){
      replicated.data <- data.frame(summed.data[1:cols])
      next
    }
    data <- data.frame(
        summed.data[seq((iii*cols)-(cols-1),(iii*cols))]
        )

    a = as.character(replicated.data[,length(replicated.data)])
    b = as.character(data$Key)
    matches <- match(a, b)
    matches <- matches[!is.na(matches)]
    data <- data[matches,]

    a = as.character(replicated.data[,length(replicated.data) ])
    b = as.character(data$Key)
    matches <- match(b, a)
    matches <- matches[!is.na(matches)]
    replicated.data <- replicated.data[matches,]

    replicated.data <- data.frame(replicated.data, data)
  }
  return(replicated.data)
}

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
  row.names=T,
  qmethod = "double"
)
