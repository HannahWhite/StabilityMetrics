############################################################
### Calculating EVI anomalies at multiple spatial scales ###
############################################################

### Hannah White 18.02.2019

## Uses MODIS data processed by Jon Yearsley and uploaded to dropbox on 18.02.2019
## hit memory errors when loading 500m data

rm(list=ls())



### 500 m x 500 m

## Load data

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_500m_all.RData')


## Calculate anomaly
# Create 18 year average baseline for  EVI


month.list = sort(unique(month))    

# For each month calculate the 18 year average and then the anomaly
evi.baseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
evi.sdbaseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
for (i in 1:12) {
  evi.baseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
  evi.sdbaseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
}


# Calculate an anomaly (deviation from baseline divided by standard deviation for each 250m square)
dm = dim(evi.median)

evi.anomaly = array(NA, dim=dm)


for (i in 1:dim(evi.median)[3]) {
  evi.anomaly[,,i] = (evi.median[,,i] - evi.baseline[,,as.numeric(month[i])])/evi.sdbaseline[,,as.numeric(month[i])]
}

save(evi.anomaly, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly500.RData')




### 1 km x 1 km


## Load data

rm(list=ls())

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')


## Calculate anomaly
# Create 18 year average baseline for  EVI


month.list = sort(unique(month))    

# For each month calculate the 18 year average and then the anomaly
evi.baseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
evi.sdbaseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
for (i in 1:12) {
  evi.baseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
  evi.sdbaseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
}


# Calculate an anomaly (deviation from baseline divided by standard deviation for each 250m square)
dm = dim(evi.median)


evi.anomaly = array(NA, dim=dm)


for (i in 1:dim(evi.median)[3]) {
  evi.anomaly[,,i] = (evi.median[,,i] - evi.baseline[,,as.numeric(month[i])])/evi.sdbaseline[,,as.numeric(month[i])]
}

save(evi.anomaly, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly1km.RData')

rm(list=ls())



### 10 km x 10 km


load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_10km_all.RData')


## Calculate anomaly
# Create 18 year average baseline for  EVI


month.list = sort(unique(month))    

# For each month calculate the 18 year average and then the anomaly
evi.baseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
evi.sdbaseline = array(NA, dim=c(dim(evi.median)[1:2], 12))
for (i in 1:12) {
  evi.baseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=mean, na.rm=TRUE)
  evi.sdbaseline[,,i] =   apply(evi.median[,,month%in%month.list[i]], MARGIN=c(1,2), FUN=sd, na.rm=TRUE)
}


# Calculate an anomaly (deviation from baseline divided by standard deviation for each 250m square)
dm = dim(evi.median)


evi.anomaly = array(NA, dim=dm)


for (i in 1:dim(evi.median)[3]) {
  evi.anomaly[,,i] = (evi.median[,,i] - evi.baseline[,,as.numeric(month[i])])/evi.sdbaseline[,,as.numeric(month[i])]
}

save(evi.anomaly, file = 'J:\\Postdoc Grassland Resilience\\MODIS6\\Anomalies_18.02.2019\\evi.anomaly10km.RData')

rm(list=ls())









