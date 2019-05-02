#####################################################################################################
### Calculating resilience metrics and timings using events > 2 standard deviations at 1 km Scale ###
#####################################################################################################

### Hannah White 18.04.2019 
### Edited 25.04.2019

rm(list=ls())

### Load data

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')
load('J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly1km.RData')

load('J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly1km.cent.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


modis.dates <- as.Date(modis.dates)

################ Resilience metric calculations #################

## index
ind2 <- year>=2003 & year<=2019

##### Resistance


julian.date2 <- julian.date[ind2]
modis.dates2 <- modis.dates[ind2]

anomaly.03 <- evi.anomaly[,,ind2] # evi anomaly array just for years 2003-2019

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
      ind.sd <- which(anomaly.03[i, j, ]<=-2) # finds anomalies which deviate more than one standard deviation from the zero baseline
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

evi.magnitude.median <- apply(evi.all.events, c(1,2), median, na.rm = TRUE) # calculates median of all drops below -2 sd

save(evi.magnitude, evi.magnitude.median, evi.events, evi.dates, evi.all.events, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km2sd.RData')

evi.dates2sd <- evi.dates
save(evi.dates2sd, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\evi.dates2sd.RData')

## Calculate and evi.drop which is the biggest drop in each site

evi.drop <- apply(evi.all.events, c(1,2), min, na.rm = TRUE) # has to be min as currently negative


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Still need to do stuff here #

## Firstly for variability and mean measures across time periods for others

# Get measures ainto dataframe


## resistance
evi.mag <- evi.magnitude[which(is.finite(evi.magnitude))]
eastings.mag <- eastings[which(is.finite(evi.magnitude))]
northings.mag <- northings[which(is.finite(evi.magnitude))]

evi.mag <- data.frame(evi.mag, eastings = eastings.mag, northings = northings.mag)

## median resistance
evi.magmed <- evi.magnitude.median[which(is.finite(evi.magnitude.median))]
eastings.magmed <- eastings[which(is.finite(evi.magnitude.median))]
northings.magmed <- northings[which(is.finite(evi.magnitude.median))]

evi.magmed <- data.frame(evi.magmed = evi.magmed, eastings = eastings.magmed, northings = northings.magmed)

## max drop
evi.maxdrop <- evi.drop[which(is.finite(evi.drop))]
eastings.maxdrop <- eastings[which(is.finite(evi.drop))]
northings.maxdrop <- northings[which(is.finite(evi.drop))]

evi.maxdrop <- data.frame(evi.drop = evi.maxdrop, eastings = eastings.maxdrop, northings = northings.maxdrop)





# change magnitude to inverse of the absolute number

stability1km$evi.mag <- 1/abs(stability1km$evi.mag)
stability1km$evi.magmed <- 1/abs(stability1km$evi.magmed)
stability1km$evi.drop <- 1/abs(stability1km$evi.drop)

write.csv(stability1km, 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\stability1km.csv', row.names = FALSE)



