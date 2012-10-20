# this file should be source()-able to import this function
library(caTools)

# data = c(1,32,33,4,3,4,3,42,2,3,4,34,34344,2,3,4,32,4,3,4,4,4,43,343,2,3,8,98,8,7,3)
# plot (data)#,ylim=c(-50000,50000))
# lines(1:length(data), kmean)
# lines(1:length(data), upper)
# lines(1:length(data), lower)
# points(data_nopeaks, pch="+")
# data <- data_nopeaks
# data[data == max(data)] = NA 
# plot(1:length(data),km.runsd(data))
# plot(1:length(data),runsd(data,5))

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

# plot(1:length(data), data)
# peak.pick(data)
# plot(1:length(data), peak.convert(data))
# peak.convert(data)
# plot(1:length(data), peaks.only(data))
# peaks.only(data)