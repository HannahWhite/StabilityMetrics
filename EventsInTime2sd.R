########################################################
### Investigating change in 2 sigma events over time ###
########################################################

### Hannah White 17.06.2019


load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


#load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km.RData')

## read in terrestrial coordinates

## read in location of terrestrial hectads

corine1km.terr <- read.csv('J:/Postdoc Grassland Resilience/LandCoverData/corine1km.terr.csv', header = TRUE)

coords <- corine1km.terr[,1:2]

load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\evi.dates2sd.RData')

### ggplot barplot
library(ggplot2)


#### All events

## stack values
dd <- evi.dates2sd[which(!is.na(evi.dates2sd))]

# extract just terrestrial sites

east.array <- array(data = NA, dim = dim(evi.dates2sd)) # create empty array to match dimensions of evi.dates so that then can use this to subset
for (i in 1:158){
  east.array[,,i] <- eastings
}
ee <- east.array[which(!is.na(evi.dates2sd))]

north.array <- array(data = NA, dim = dim(evi.dates2sd))
for (i in 1:158){
  north.array[,,i] <- northings
}
nn <- north.array[which(!is.na(evi.dates2sd))]


dd <- data.frame(dd, ee, nn)
dd$ee.nn <- paste(dd$ee, dd$nn, sep = '_')

coords$ee.nn <- paste(coords$east, coords$north, sep = '_')

dd <- dd[(which(dd$ee.nn %in% coords$ee.nn)),] # extracts rows of dd where eastings and northings are in coords df which is just terrestrial squares


## calculations

tab <- table(cut(dd$dd, 'year'))

# take out 2019
tab <- tab[-17]


year.df <- data.frame(year=format(as.Date(names(tab)), '%Y'), frequency=as.vector(tab))


m1 <- lm(frequency~as.numeric(year), data = year.df)
summary(m1)

m2 <- glm(frequency~as.numeric(year), data = year.df, family ='poisson')
summary(m2) # massively overdispersed


library(MASS)
m3 <- glm.nb(frequency~as.numeric(year), data = year.df, link = 'log')
summary(m3)

plot(residuals(m3, type = 'deviance') ~ fitted(m3))

m4 <- glm.nb(frequency~1, data = year.df, link = 'log') # 69438.03

## ggplot2

p <- ggplot(year.df, aes(x = year, y = frequency)) +  geom_point(fill = '#00AFBB', colour = 'black', size = 4, shape = 21, stroke = 1.5)   #geom_bar(stat='identity')
p <- p + xlab('Year') + ylab('Frequency') + ylim(60000, 90950)
p <- p + geom_hline(yintercept = 83403, colour = 'red', linetype = 'dashed')
p <- p + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
               panel.grid = element_blank(), panel.background = element_blank(),
               axis.line = element_line(colour = 'black'))

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\EventsInTime1km.tiff', width = 14, height = 14, units = 'in', res = 360)
p
dev.off()


## base plotting with plotrix to get axis break
library(plotrix)

par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(frequency ~ as.Date(year, '%Y'), data = year.df, pch = 16, col = 'black',
     ylim = c(50000, 80000), cex = 3, xlab = 'Year', ylab = 'Number of two sigma events', cex.axis = 2, cex.lab =2)
abline(a = 69438, b = 0, col = 'red', lty = 2, lwd = 3)

axis.break(axis = 2, breakpos = 49500, style = 'gap', brw = 0.015)

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\Sigma2Events\\EventsInTime1km_base.tiff', width = 14, height = 14, units = 'in', 
     res = 360, compression = 'lzw')
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(frequency ~ as.Date(year, '%Y'), data = year.df, pch = 16, col = 'black',
     ylim = c(50000, 80000), cex = 3, xlab = 'Year', ylab = 'Number of two sigma events', cex.axis = 2, cex.lab =2)
abline(a = 69438, b = 0, col = 'red', lty = 2, lwd = 3)

axis.break(axis = 2, breakpos = 49500, style = 'gap', brw = 0.015)
dev.off()

