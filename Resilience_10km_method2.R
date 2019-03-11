#######################################################################################
### Calculating resilience metrics using biggest drop window methods at 10 km Scale ###
#######################################################################################

### Hannah White 14.02.2019
### Edited 18.02.2019 to calculate return times 


### This is a development of previous code where now, large perturbations in the anomalies are found and 180 days each side are excluded

rm(list=ls())



## read in anomaly data

load('J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly10km.RData')


## read in EVI data so that have Julian dates

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_10km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


modis.dates <- as.Date(modis.dates)

################ Resilience metric calculations #################

## index
ind2 <- year>=2003 & year<=2019

##### Variability

evi.variability <- array(NA, dim =c(dim(eastings)))


## calculate variability from 2000-2019

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

evi.all.events <- array(NA, dim=c(dim(eastings), 158)) # made 3rd dimension 158 as maximum number of resistance events is 157

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


### quickly plot some histograms

### ggplot histogram
#library(ggplot2)
#library(scales)
#library(lubridate)


#### All events

## stack values
#dd <- evi.dates[which(!is.na(evi.dates))]
#ddnum <- as.numeric(dd)

#dd.df <- data.frame(dd, ddnum)

#bin = 30

#p <- ggplot(dd.df, aes(dd, ..count..))
#p <- p + geom_histogram(data = dd.df, binwidth = bin)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
#p <- p + scale_x_date(breaks = waiver(),
#                      date_breaks = '1 year',
#                      labels = date_format("%Y-%m"),
#                      limits = c(as.Date("2003-01-01"), 
#                                 as.Date("2019-01-25")))
#p <- p + xlab('')
#p <- p + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
#p <- p + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
#               panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

#p


#### Largest events

#evi.dates.big <- evi.dates[,,1]
#ddbig <- evi.dates.big[which(!is.na(evi.dates.big))]
#ddnumbig <- as.numeric(ddbig)

#ddbig.df <- data.frame(ddbig, ddnumbig)


#p2 <- ggplot(ddbig.df, aes(x = ddbig, y = ..count..))
#p2 <- p2 + geom_histogram(data = ddbig.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
#p2 <- p2 + scale_x_date(breaks = waiver(),
#                      date_breaks = '1 year',
#                      labels = date_format("%Y-%m"),
#                      limits = c(as.Date("2003-01-01"), 
#                                 as.Date("2019-01-25")))
#p2 <- p2 + xlab('')
#p2 <- p2 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
#p2 <- p2 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
#                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

#p2

#library(cowplot)

#tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\events10km.tiff', width = 14, height = 7, units = 'in', res = 360)
#plot_grid(p, p2, nrow = 1, labels = c('a)', 'b)'))
#dev.off()


### Check for any seasonality using decompose

date.count <- table(ddbig.df[,1])

date.count <- data.frame(date.count)
names(date.count) <- c('date', 'count')

decompose(date.count, type = ("additive"), filter = NULL)

data.ts <- ts(date.count[,2], start = c(2003,1), frequency = 12) # frequency isn't quite right date wise but get the gist
plot(decompose(data.ts, type = 'additive'))


##### Recovery time


### First need to calculate moving windows

# Create empty array to hold moving window averages where window size is 6 units

mov.avg <- array(NA, dim = c(dim(eastings),(dim(anomaly.03)[3])-5))


for (kk in 1:734){ #begins moving window at each separate time point
  window.frame <- anomaly.03[,,kk:(kk+5)] # moving window of 6 units
  mov.avg[,,kk] <- apply(window.frame, c(1,2), mean, na.rm = TRUE)
  
}


# Return times  and rates - calculated simultaneously

evi.all.rectime <- array(NA, dim=c(dim(eastings), 158))


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
            
            if (start<=733) {
          
           return.time <- which(mov.avg[ii, jj, (start+1):734]>=0)[1] # calculates recovery time in units of time windows
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


##### Save all data in same place

#save(evi.all.events, evi.dates, evi.dates.big, evi.events, evi.magnitude, evi.all.rectime, evi.all.rate, evi.rectime, evi.rate, evi.variability,
#     file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_10km.RData')


### Cross Correlations

## Calculate and evi.drop which is the biggest drop in each site

evi.drop <- apply(evi.all.events, c(1,2), min, na.rm = TRUE) # has to be min as currently negative



## Firstly for variability and mean measures across time periods for others

# Get mean measures and variability into dataframe

## var
evi.var <- evi.variability[which(is.finite(evi.variability))]
eastings.var <- eastings[which(is.finite(evi.variability))]
northings.var <- northings[which(is.finite(evi.variability))]

evi.var <- data.frame(evi.var, eastings = eastings.var, northings = northings.var)


## resistance
evi.mag <- evi.magnitude[which(is.finite(evi.magnitude))]
eastings.mag <- eastings[which(is.finite(evi.magnitude))]
northings.mag <- northings[which(is.finite(evi.magnitude))]

evi.mag <- data.frame(evi.mag, eastings = eastings.mag, northings = northings.mag)

## max drop
evi.maxdrop <- evi.drop[which(is.finite(evi.drop))]
eastings.maxdrop <- eastings[which(is.finite(evi.drop))]
northings.maxdrop <- northings[which(is.finite(evi.drop))]

evi.maxdrop <- data.frame(evi.drop = evi.maxdrop, eastings = eastings.maxdrop, northings = northings.maxdrop)

## recovery time
evi.rec <- evi.rectime[which(is.finite(evi.rectime))]
eastings.rec <- eastings[which(is.finite(evi.rectime))]
northings.rec <- northings[which(is.finite(evi.rectime))]

evi.rec <- data.frame(evi.rec, eastings = eastings.rec, northings = northings.rec)


## rate
evi.recrate <- evi.rate[which(is.finite(evi.rate))]
eastings.recrate <- eastings[which(is.finite(evi.rate))]
northings.recrate <- northings[which(is.finite(evi.rate))]

evi.recrate <- data.frame(evi.recrate, eastings = eastings.recrate, northings = northings.recrate)


### Merge

stability10km <- merge(evi.var, evi.mag, by = c('eastings', 'northings'), all = TRUE)
stability10km <- merge(stability10km, evi.rec, by = c('eastings', 'northings'), all = TRUE)
stability10km <- merge(stability10km, evi.recrate, by = c('eastings', 'northings'), all = TRUE)
stability10km <- merge(stability10km, evi.maxdrop, by = c('eastings', 'northings'), all = TRUE)


# change magnitude to inverse of the absolute number

stability10km$evi.mag <- 1/abs(stability10km$evi.mag)
stability10km$evi.drop <- 1/abs(stability10km$evi.drop)

cor(stability10km[,3:7], use = 'na.or.complete', method = 'spearman') # use Spearman's as many correlations appear non linear (see pair graph)


#panel.cor <- function(x, y, ...)
#{
#  par(usr = c(0, 1, 0, 1))
#  txt <- as.character(format(cor(x, y, use = 'na.or.complete', method = 'spearman'), digits=2))
#  text(0.5, 0.5, txt, cex = 4)
#}

#tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\correlations10km.tiff', height = 1000, width = 1000)
#pairs(stability10km[,3:7], labels = c('variability', 'resistance', 'recovery time', 'recovery rate', 'biggest perturbation'), 
#      upper.panel = panel.cor) ## need to change colour of some of them as non-sig
#dev.off()

library(Hmisc)

rcorr(cbind(stability10km[,3], stability10km[,4], stability10km[,5], stability10km[,6], stability10km[,7]),  type = 'spearman') ## all significant







