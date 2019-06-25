#####################################################################################################
### Analyses and graphics for stability metrics at 1 km resolution using sigma 2 events and means ###
#####################################################################################################

### Hannah White 20.06.2019

library(ggplot2)
library(cowplot)
library(raster)
library(viridis)
library(scales)


## load stability measures

stability1km <- read.csv('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\stability1km2sd.csv', header = TRUE)

## read in location of terrestrial hectads

corine1km.terr <- read.csv('J:/Postdoc Grassland Resilience/LandCoverData/corine1km.terr.csv', header = TRUE)

coords <- corine1km.terr[,1:2]


stability1km <- merge(stability1km, coords, by.x = c('eastings', 'northings'), by.y = c('east', 'north'), all.x = FALSE, all.y = TRUE)

## Change recovery time to days rather than units

stability1km$evi.recmed <- stability1km$evi.recmed*8
stability1km$evi.rec <- stability1km$evi.rec*8


stability.red <- stability1km[, c(1:2, 12, 3, 5, 6, 13:14)] # variability, mag, rec, recrate, long, slow
#stability.red.log <- data.frame(stability.red[,1:2], log(stability.red[,3:8]))

stability.redless <- stability.red

# remove data points which will change scale on maps

stability.redless$evi.var[which(stability.redless$evi.var>0.02)] <- NA #changes 140
stability.redless$evi.mag[which(stability.redless$evi.mag>0.65)] <- NA #changes 194
stability.redless$evi.rec[which(stability.redless$evi.rec>100)] <- NA #changes 145
stability.redless$evi.long[which(stability.redless$evi.long>100)] <- NA #changes 246
stability.redless$evi.slow[which(stability.redless$evi.slow>0.002)] <- NA #changes 217

### Maps

dfr <- rasterFromXYZ(stability.redless)  #Convert first two columns as lon-lat and third as value    


tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\Sigma2Events\\maps1kmMEAN.tiff', 
     height = 7, width = 14, units = 'in', res = 360, compression = 'lzw') #may need to adjust cut off point
par(oma = c(0,0,0,6))
plot(dfr, box = FALSE, axes = FALSE, 
     main = c('Variability', 'Resistance', 'Recovery Time', 'Recovery Rate', 'Longest Return', 'Slowest Return'),
     cex.main = 2, legend.width = 2, 
     axis.args=list(cex.axis=2.5))
dev.off()




### Cross-correlations


#pairs(stability.red[,3:8])


cor(stability.red[,3:8], use = 'na.or.complete', method = 'spearman') # use Spearman's as many correlations appear non linear (see pair graph)

#pairs(stability.red[,3:8]) # low correlations, so seem to capturing different aspects of stability

#panel.cor <- function(x, y, ...)
#{
#  par(usr = c(0, 1, 0, 1))
#  txt <- as.character(format(cor(x, y, use = 'na.or.complete', method = 'spearman'), digits=2))
#  text(0.5, 0.5, txt, cex = 4)
#}

#tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\correlations1km.tiff', height = 1000, width = 1000)
#pairs(stability.red[,3:8], labels = c('variability', 'resistance', 'recovery time', 'recovery rate', 'longest return', 'slowest return'), 
#      upper.panel = panel.cor, cex.axis = 2) ## need to change colour of some of them as non-sig
#dev.off()

library(Hmisc)

rcorr(cbind(stability.red[,3], stability.red[,4], stability.red[,5], stability.red[,6], stability.red[,7], stability.red[,8]),  type = 'spearman') ## all significant



# pairs plot with histogram
library(psych)

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\Sigma2Events\\correlations1kmhistMEAN.tiff', 
     height = 14, width = 14, units = 'in', res = 360, compression = 'lzw')
pairs.panels(stability.red[,3:8],
             labels = c('variability', 'resistance', 'recovery time', 'recovery rate', 'longest return', 'slowest return'),
             method = "spearman", # correlation method
             smooth = FALSE, # don't draw loess smooth
             hist.col = "#00AFBB",
             rug = FALSE,
             density = TRUE, # show density plots
             ellipses = FALSE,
             cex.axis = 2
)
dev.off()

### Timing of events


load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km.RData')

load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\evi.dates2sd.RData')

### ggplot histogram
library(ggplot2)
library(scales)
library(lubridate)

#### All events

## stack values
dd <- evi.dates[which(!is.na(evi.dates))]
ddnum <- as.numeric(dd)

dd.df <- data.frame(dd, ddnum)


## eastings and northings
deast <- eastings[which(!is.na(evi.dates[,,1]))]

for (i in 2:dim(evi.dates)[3]){
  east.next <- eastings[which(!is.na(evi.dates[,,i]))]
  deast <- c(deast, east.next)
}


dnorth <- northings[which(!is.na(evi.dates[,,1]))]

for (j in 2:dim(evi.dates)[3]){
  north.next <- northings[which(!is.na(evi.dates[,,j]))]
  dnorth <- c(dnorth, north.next)
}

dd.df$eastings <- deast
dd.df$northings <- dnorth

coords$comb <- paste(coords$east, coords$north, sep = '.')

dd.df$comb <- paste(dd.df$eastings, dd.df$northings, sep = '.')


dd.df <- merge(dd.df, coords, by = 'comb', all.x = FALSE, all.y = TRUE) # extracts those with > 50% land

dd.df$eventsd <- 'one'


## Events >2sd

dd2sd <- evi.dates2sd[which(!is.na(evi.dates2sd))]
ddnum2sd <- as.numeric(dd2sd)

dd2sd.df <- data.frame(dd2sd, ddnum2sd)

deast2 <- eastings[which(!is.na(evi.dates2sd[,,1]))]

for (i in 2:dim(evi.dates2sd)[3]){
  east.next2sd <- eastings[which(!is.na(evi.dates2sd[,,i]))]
  deast2 <- c(deast2, east.next2sd)
}

dnorth2 <- northings[which(!is.na(evi.dates2sd[,,1]))]

for (j in 2:dim(evi.dates2sd)[3]){
  north.next2sd <- northings[which(!is.na(evi.dates2sd[,,j]))]
  dnorth2 <- c(dnorth2, north.next2sd)
}

dd2sd.df$eastings <- deast2
dd2sd.df$northings <- dnorth2

dd2sd.df$comb <- paste(dd2sd.df$eastings, dd2sd.df$northings, sep = '.')

#mm <- merge(dd2sd.df, coords, by = 'comb', all.x = FALSE, all.y = FALSE) #use this line to calculate % over 2sd

dd2sd.df <- merge(dd2sd.df, coords, by = 'comb', all.x = FALSE, all.y = TRUE)


dd2sd.df$eventsd <- 'two'


names(dd2sd.df)[2:3] <- c('dd', 'ddnum')


dd.all <- rbind(dd.df, dd2sd.df)


### Plot



bin = 30

p <- ggplot(dd.df, aes(dd, fill = eventsd))
p <- p + geom_histogram(data = dd.all, binwidth = bin, aes(y = ..count..), position = 'identity')

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p <- p + scale_x_date(breaks = waiver(),
                      date_breaks = '1 year',
                      labels = date_format("%Y-%m"),
                      limits = c(as.Date("2003-01-01"), 
                                 as.Date("2019-01-25")))
p <- p + xlab('')
p <- p + scale_fill_manual(values = c('#00AFBB', 'grey20'))

# from here, format at ease
p <- p + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 15), axis.text.y = element_text(size = 14), legend.position ='none',
               panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p






#### Largest events

evi.dates.big <- evi.dates[,,1]
ddbig <- evi.dates.big[which(!is.na(evi.dates.big))]
ddnumbig <- as.numeric(ddbig)

ddbig.df <- data.frame(ddbig, ddnumbig)

deast <- eastings[which(!is.na(evi.dates[,,1]))]
dnorth <- northings[which(!is.na(evi.dates[,,1]))]

ddbig.df$east <- deast
ddbig.df$north <- dnorth

ddbig.df <- merge(ddbig.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE) # extracts those with >50% land

#ddbig.df$eventsd <- 'one'

## 2 sd events

#evi.dates.big2sd <- evi.dates2sd[,,1]
#ddbig2sd <- evi.dates.big2sd[which(!is.na(evi.dates.big2sd))]
#ddnumbig2sd <- as.numeric(ddbig2sd)

#ddbig2sd.df <- data.frame(ddbig2sd, ddnumbig2sd)

#deast2 <- eastings[which(!is.na(evi.dates.big2sd))]
#dnorth2 <- northings[which(!is.na(evi.dates.big2sd))]

#ddbig2sd.df$east <- deast2
#ddbig2sd.df$north <- dnorth2

#ddbig2sd.df <- merge(ddbig2sd.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE)

#ddbig2sd.df$eventsd <- 'two'

#names(ddbig2sd.df)[3:4] <- c('ddbig', 'ddnumbig')


#ddbig.all <- rbind(ddbig.df, ddbig2sd.df)


## plot

p2 <- ggplot(ddbig.df, aes(x = ddbig))
p2 <- p2 + geom_histogram(data = ddbig.df, binwidth = 30, aes(y = ..count..), position = 'identity')

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p2 <- p2 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p2 <- p2 + xlab('')
#p2 <- p2 + scale_fill_manual(values = c('grey70', 'grey20'))

# from here, format at ease
p2 <- p2 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 15), axis.text.y = element_text(size = 14), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p2

library(cowplot)

#tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\events10km.tiff', width = 14, height = 7, units = 'in', res = 360)
#plot_grid(p, p2, nrow = 1, labels = c('a)', 'b)'))
#dev.off()


### Longest return rates



ddlong <- evi.longest.date[which(!is.na(evi.longest.date))]
ddnumlong <- as.numeric(ddlong)

ddlong.df <- data.frame(ddlong, ddnumlong)

deast <- eastings[which(!is.na(evi.longest.date))]
dnorth <- northings[which(!is.na(evi.longest.date))]

ddlong.df$east <- deast
ddlong.df$north <- dnorth

ddlong.df <- merge(ddlong.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE) # extracts those with >50% land



p3 <- ggplot(ddlong.df, aes(x = ddlong, y = ..count..))
p3 <- p3 + geom_histogram(data = ddlong.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p3 <- p3 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p3 <- p3 + xlab('')
#p3 <- p3 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p3 <- p3 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 15), axis.text.y = element_text(size = 14), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p3

## Slowest events

ddslow <- evi.slowest.date[which(!is.na(evi.slowest.date))]
ddnumslow <- as.numeric(ddslow)

ddslow.df <- data.frame(ddslow, ddnumslow)

deast <- eastings[which(!is.na(evi.slowest.date))]
dnorth <- northings[which(!is.na(evi.slowest.date))]

ddslow.df$east <- deast
ddslow.df$north <- dnorth

ddslow.df <- merge(ddslow.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE)



p4 <- ggplot(ddslow.df, aes(x = ddslow, y = ..count..))
p4 <- p4 + geom_histogram(data = ddslow.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p4 <- p4 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p4 <- p4 + xlab('')
#p4 <- p4 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p4 <- p4 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 15), axis.text.y = element_text(size = 14), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p4



tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\events1kmextra.tiff', width = 14, height = 14, units = 'in', res = 360, compression = 'lzw')
plot_grid(p, p2, p3, p4, nrow = 2, ncol = 2, labels = c('a)', 'b)', 'c)', 'd)'))
dev.off()




