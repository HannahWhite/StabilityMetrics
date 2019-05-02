###########################################################################
### Calculating slowest recovery rate and longest recovery time at 1 km ###
###########################################################################

### Hannah White 28.03.2019
### Edited 02.04.2019 to work out when these are occurring

# Calculates the longest recovery time and slowest return rate of all of the EVI anomalies.

rm(list=ls())



## read in anomaly data

load('J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly1km.RData')


## read in EVI data so that have Julian dates

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


modis.dates <- as.Date(modis.dates)

################ Resilience metric calculations #################

## index
ind2 <- year>=2003 & year<=2019


julian.date2 <- julian.date[ind2]
modis.dates2 <- modis.dates[ind2]

anomaly.03 <- evi.anomaly[,,ind2] # evi anomaly array just for years 2003-2017


### First need to calculate moving windows

# Create empty array to hold moving window averages where window size is 6 units

mov.avg <- array(NA, dim = c(dim(eastings),(dim(anomaly.03)[3])-5))


for (kk in 1:734){ #begins moving window at each separate time point
  window.frame <- anomaly.03[,,kk:(kk+5)] # moving window of 6 units
  mov.avg[,,kk] <- apply(window.frame, c(1,2), mean, na.rm = TRUE)
  
}


# Return times  and rates - calculated simultaneously

recovery.full <- array(NA, dim=c(dim(eastings), (dim(anomaly.03)[3])-5))


rate.full <- array(NA, dim=c(dim(eastings), (dim(anomaly.03)[3])-5))


### need to match dates

for (ii in 1:dim(eastings)[1]){
  for (jj in 1:dim(eastings)[2]){
    if (any(!is.na(anomaly.03[ii,jj,]))){
      
      #event.length <- length(which(is.finite(evi.all.events[ii,jj,]))) ## number of 'resistance' events
      #date.event <- evi.dates[ii, jj, ] # need to extract date so that know when to start searching moving windows
      #date.event.finite <- date.event[which(!is.na(date.event))] 
      
      
      for (mm in 1:(dim(anomaly.03)[3]-6)){ # need to subtract 6 so that last one has dates to search within
        
        
        
        return.time <- which(mov.avg[ii, jj, (mm+1):734]>=0)[1] # calculates recovery time in units of time windows
        # needs to be start + 1 so that event window isn't counted
        recovery.full[ii, jj, mm] <- return.time
        
        rate.full[ii, jj, mm] <- abs(anomaly.03[ii, jj, mm])/return.time # uses absolute so that rates are positive 
        
        
      } 
      
      
    }
    
  }
  
}


evi.longest.rec <- apply(recovery.full, c(1,2), max, na.rm = TRUE)

evi.slowest.rate <- apply(rate.full, c(1,2), min, na.rm = TRUE)



which.max.na <- function(x){
  ifelse(length(which(is.na(x)))==734, NA, which.max(x))
}


which.min.na <- function(x){
  ifelse(length(which(is.na(x)))==734, NA, which.min(x))
}



evi.longest.which <- apply(recovery.full, c(1,2), function(x) which.max.na(x))

evi.slowest.which <- apply(rate.full, c(1,2), function(x) which.min.na(x))


evi.longest.date <- array(NA, dim = c(dim(evi.longest.which)))
class(evi.longest.date) <- 'Date'

evi.slowest.date <- array(NA, dim = c(dim(evi.slowest.which)))
class(evi.slowest.date) <- 'Date'

for(i in 1:dim(eastings)[1]){
  for(j in 1:dim(eastings)[2]){
    
    ind.long <- evi.longest.which[i,j] 
    evi.longest.date[i,j] <- modis.dates2[ind.long]
    
    
    ind.slow <- evi.slowest.which[i,j]
    evi.slowest.date[i,j] <- modis.dates2[ind.slow]
    
  }
}



