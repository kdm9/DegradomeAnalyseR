library(plyr)
# setwd("/home/kevin/UniWork/Pogson Lab/Degradome/PetePARE/paresnip_results")

p.list <- list(
    list(
        c("paresnip_1.tab","WT"),
        c("paresnip_2.tab","WT"),
        c("paresnip_3.tab","WT")
        ),
    list(
        c("paresnip_4.tab","xrn4"),
        c("paresnip_5.tab","xrn4"),
        c("paresnip_6.tab","xrn4")
        ),
    list(
        c("paresnip_7.tab","alx8"),
        c("paresnip_8.tab","alx8"),
        c("paresnip_9.tab","alx8")
        )
    )
params = list(
      c("paresnip_1.tab","WT"),
      c("paresnip_2.tab","WT"),
      c("paresnip_3.tab","WT")
    )


cross.validate <- function (params){
  counter = 0
  summed.data = list()
  replicated.data = list()
  for (param in params){
    counter = counter + 1
    data <- read.delim(param[1],header=T, quote = '"', sep='\t')
    agi.list <- unlist(strsplit(as.character(data$Gene), "|", fixed=T))
    data$AGI <- matrix(agi.list[grep("^AT.G", agi.list)])

    data.summarised <- ddply(
      data[data$Category<3,],
      .(AGI, Cleavage.Position, Category),
      summarize,
      FA=max(Fragment.Abundance),
      SRA=max(Short.Read.Abundance),
      Duplex=names(which.max(table(as.character(Duplex)))) #most common Duplex
    )

    data.summarised$Key <- paste(
      data.summarised$AGI,
      data.summarised$Cleavage.Position,
      data.summarised$Category
      )

    data.summarised <- data.summarised[with(data.summarised, order(Key)), ]
    names(data.summarised) <- paste(names(data.summarised), counter)
    summed.data <- c(summed.data, data.summarised)
    print(length(summed.data))

  }
  print(names(summed.data))
  for (iii in c(1:(counter))){
    cols = 7
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


for (set in p.list){
    fn = paste(set[[1]][2], "conserved.tab", sep="_")
    print(fn)
    conserved <- data.frame(cross.validate(set))

    write.table(
        conserved,
        file=fn,
        sep="\t",
        col.names=T,
        row.names=T,
        qmethod = "double"
        )
}
