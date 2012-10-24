library(rgl)
args <- commandArgs(trailingOnly=T)

paresnip.summarised <- read.csv(args[1])

pdf(paste(args[2], "corelations.pdf", sep="_"))

# FA vs SRA relationship (there is none)
plot(paresnip.summarised$FA.mean, paresnip.summarised$SRA.mean, main="SRA vs FA")
plot(paresnip.summarised$FA.mean, paresnip.summarised$SRA.mean, ylim=c(0,5000),xlim=c(0,5000), main="SRA vs FA clipped")
FAvSRA <- lm(paresnip.summarised$FA.mean ~ paresnip.summarised$SRA.mean)
summary(FAvSRA)
abline(FAvSRA)

# transforms
plot(1/(paresnip.summarised$FA.mean), paresnip.summarised$SRA.mean, main="SRA vs 1/FA")
inv_FAvSRA <- lm(1/paresnip.summarised$FA.mean ~ paresnip.summarised$SRA.mean)
summary(inv_FAvSRA)
abline(inv_FAvSRA)

plot(log(paresnip.summarised$FA.mean), paresnip.summarised$SRA.mean, main="SRA vs ln(FA)")
log_FAvSRA <- lm(log(paresnip.summarised$FA.mean) ~ paresnip.summarised$SRA.mean)
summary(log_FAvSRA)
abline(log_FAvSRA)

# Look at inter-sample replication in FA and SRA
plot((paresnip.summarised$FA.1), (paresnip.summarised$FA.2), main="FA 1v2")
abline(a=0,b=1)
plot((paresnip.summarised$FA.3), (paresnip.summarised$FA.2), main="FA 3v2")
abline(a=0,b=1)
plot((paresnip.summarised$FA.1), (paresnip.summarised$FA.3), main="FA 1v3")
abline(a=0,b=1)

plot((paresnip.summarised$SRA.1), (paresnip.summarised$SRA.2), main="SRA 1v2")
abline(a=0,b=1)
plot((paresnip.summarised$SRA.3), (paresnip.summarised$SRA.2), main="SRA 3v2")
abline(a=0,b=1)
plot((paresnip.summarised$SRA.1), (paresnip.summarised$SRA.3), main="SRA 1v3")
abline(a=0,b=1)

dev.off()

#plot3d not working on desktop, but the below should work
plot3d(paresnip.summarised$FA.1, paresnip.summarised$FA.2, paresnip.summarised$FA.3, main="FA between samples")
plot3d(paresnip.summarised$SRA.1 ,paresnip.summarised$SRA.2, paresnip.summarised$SRA.3, main="SRA between samples")
