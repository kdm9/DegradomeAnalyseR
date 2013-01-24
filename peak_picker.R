# DegradomeAnalyseR is copyright 2012 Kevin Murray, and is licensed under the
# GPLv3 License

# this file should be source()-able to import this function
library(caTools)

km.runsd <- function(data){
  data.runsd <- numeric(length=length(data))
  for (i in 1:length(data)){
    data.na <- data
    data.na[i] = NA
    data.na.runsd <- runsd(data.na,5)
    data.runsd[i] = data.na.runsd[i]
  }
  return(data.runsd)
}

peak.pick <- function(data){
  kmean <- runmean(data,15)
  ksd = 1 * sd(data)
  upper <- kmean + ksd
  lower <- kmean - ksd
  peaks <- (upper < data) | (data < lower)
  return(peaks)
}

peak.convert <- function (data){
  peaks <- peak.pick(data)
  data_nopeaks <- data
  data[peaks] = NA
  kmean <- runmean(data,15)
  data_nopeaks[peaks] = kmean[peaks]
  return(data_nopeaks)
}

peaks.only <- function (data){
  peaks <- peak.pick(data)
  data_peaks <- data
  data_peaks[!peaks] = NA
  return(data_peaks)
}
