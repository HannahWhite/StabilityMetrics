######################################################################################
### Calculating resilience metrics using biggest drop window methods at 1 km Scale ###
######################################################################################

### Hannah White 18.02.2019 

rm(list=ls())

### Load data

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')
load('J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly1km.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


modis.dates <- as.Date(modis.dates)

################ Resilience metric calculations #################

## index
ind2 <- year>=2003 & year<=2019

##### Variability

evi.variability <- array(NA, dim =c(dim(eastings)))


## calculate variability from 2000-2017

for (i in 1:dim(eastings)[1]){
  for (j in 1:dim(eastings)[2]){
    if (any(is.finite(evi.anomaly[i,j,]))) evi.variability[i,j] <- var(na.omit(evi.anomaly[i,j,]))
    
  }
}


##### Resistance


julian.date2 <- julian.date[ind2]
modis.dates2 <- modis.dates[ind2]

anomaly.03 <- evi.anomaly[,,ind2] # evi anomaly array just for years 2003-2017

evi.magnitude <- array(NA, dim=c(dim(eastings)))

evi.events <- array(NA, dim=c(dim(eastings)))

evi.all.events <- array(NA, dim=c(dim(eastings), 158)) # made 3rd dimension 158 as maximum number of resistance events is 157 - actually 22 at this scale

evi.dates <- array(NA, dim=c(dim(eastings), 158))
class(evi.dates) <- 'Date' # for array to fill with dates needs to be set as Date class


### testing
#i <- 23
#j <- 18

## Function returns NA is y is within 6 months of x
magnitude.compare <- function(vec, x, y){
  return(
    ifelse(vec[y]-vec[x]<=180 & vec[y]-vec[x]>=-180, NA, order.anom[y])
  )
}

#Calculating magnitude (resistance) 

try(
  for (i in 1:dim(eastings)[1]){
    for (j in 1:dim(eastings)[2]){
      ind.sd <- which(anomaly.03[i, j, ]<=-1) # finds anomalies which deviate more than one standard deviation from the zero baseline
      anomaly.sub <- anomaly.03[i,j,ind.sd] # subsets the anomaly dataframe to just those that are greater than one standard deviation
      
      julian.ind <- julian.date2[ind.sd]  # finds the julian dates for these anomalies
      dates.ind <- modis.dates2[ind.sd]
      
      order.anom <- anomaly.sub[order(anomaly.sub)]  # orders anomalies from largest to smallest
      order.julian <- julian.ind[order(anomaly.sub)] # orders julian dates in same order as anomalies
      order.dates <- dates.ind[order(anomaly.sub)] # orders modis dates in same order as anomalies
      
      if(length(ind.sd)>0){
        close.ind <- array(NA, dim = c(length(ind.sd), length(ind.sd)))
        
        
        for (k in 1:length(ind.sd)) {
          for (l in 1:length(ind.sd)){
            
            close.ind[k, l] <- magnitude.compare(order.julian, k, l) # compares number of days between anomalies, if within 6 months each way, returns NA
          } 
        }
        
        diag(close.ind) <- order.anom # sets diagonal to anomalies rather than NA
        
        resistances <- close.ind[1,] # sets up resistances vector initially with anomalies kept that fall outside of 6 months of the largest anomaly
        
        if(length(ind.sd)>1){
          
          for (m in 2:length(ind.sd)){   ## loops through close.ind matrix diagonally and adds new NAs
            sub <- close.ind[m,]    
            na.ind <- which(is.na(sub))
            na.ind2 <- na.ind[which(na.ind > m)]
            resistances[na.ind2] <- NA
            
          } 
          
          dates <- order.dates
          dates[which(is.na(resistances))] <- NA
          
        } else {
          resistances <- resistances
          dates <- order.dates
        } 
        
        
      } else { 
        resistances <- NA
        dates <- NA
      }
      
      
      evi.all.events[i,j, ] <- c(resistances, rep(NA, 158-length(resistances))) # fills each site with all resistances plus remainder NAs
      evi.dates[i,j, ] <- c(dates, rep(NA, 158-length(dates))) 
      
      evi.events[i,j] <- ifelse(length(resistances)==length(which(is.na(resistances))), 0, length(which(!is.na(resistances)))) # if resistances is just NA then indicates there are 0 resistance events
      evi.magnitude[i,j] <- mean(resistances, na.rm = TRUE) # calculates mean of resistances for each individual pixel
    }
  }
)

save(evi.magnitude, evi.events, evi.dates, evi.all.events, evi.variability, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km.RData')


##### Recovery time


### First need to calculate moving windows

# Create empty array to hold moving window averages where window size is 6 units

mov.avg <- array(NA, dim = c(dim(eastings),(dim(anomaly.03)[3])-5))


for (kk in 1:685){ #begins moving window at each separate time point
  window.frame <- anomaly.03[,,kk:(kk+5)] # moving window of 6 units
  mov.avg[,,kk] <- apply(window.frame, c(1,2), mean, na.rm = TRUE)
  
}


# Return times  and rates - calculated simultaneously

evi.all.rectime <- array(NA, dim=c(dim(eastings), 158)) # actually could change this to 23 as max number of events at this scale is 22


evi.all.rate <- array(NA, dim=c(dim(eastings), 158))



## test
#ii <- 23
#jj <- 19

### need to match dates

for (ii in 1:dim(eastings)[1]){
  for (jj in 1:dim(eastings)[2]){
    if (any(!is.na(anomaly.03[ii,jj,]))){
      
      event.length <- length(which(is.finite(evi.all.events[ii,jj,]))) ## number of 'resistance' events
      date.event <- evi.dates[ii, jj, ] # need to extract date so that know when to start searching moving windows
      date.event.finite <- date.event[which(!is.na(date.event))] 
      
      if (event.length > 0){
        for (mm in 1:event.length){
          
          
          d.day <- date.event.finite[mm]
          start <- which(modis.dates2==d.day)
          
          if (start<=684) {
            
            return.time <- which(mov.avg[ii, jj, (start+1):685]>=0)[1] # calculates recovery time in units of time windows
            # needs to be start + 1 so that event window isn't counted
            
            evi.all.rectime[ii, jj, mm] <- return.time
            
            evi.all.rate[ii, jj, mm] <- abs(anomaly.03[ii, jj, start])/return.time # uses absolute so that rates are positive 
            
            
          } else {
            evi.all.rectime[ii, jj, mm] <- NA # gives recovery of NA if there are no events that dropped below -1 standard deviations
            evi.all.rate[ii, jj, mm] <- NA # gives recovery rate of NA if there are no events that dropped below -1 standard deviations
          }
          
          
        } 
        
        
      }
    }
    
  }
  
}


evi.rectime <- apply(evi.all.rectime, c(1,2), mean, na.rm = TRUE)  # calculates mean return time 
evi.rate <- apply(evi.all.rate, c(1,2), mean, na.rm = TRUE)

save(evi.magnitude, evi.events, evi.dates, evi.all.events, evi.variability, 
     evi.all.rectime, evi.all.rate, evi.rectime, evi.rate, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km.RData')

