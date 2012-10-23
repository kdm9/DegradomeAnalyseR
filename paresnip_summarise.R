args <- commandArgs(trailingOnly=T)
args[1] = "alx8_conserved.tab"
paresnip.replicated <- read.delim(args[1], header=T, sep='\t', quote="\"")


