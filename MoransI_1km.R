#############################################################
### Investigation of spatial autocorrelation at 1km scale ###
#############################################################

### Hannah White 05.04.2019

library(spdep)
library(spatial)
library(cowplot)

### Prepare data

## load stability measures

stability1km <- read.csv('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\stability1km.csv', header = TRUE)

## read in location of terrestrial hectads

corine1km.terr <- read.csv('J:/Postdoc Grassland Resilience/LandCoverData/corine1km.terr.csv', header = TRUE)

coords <- corine1km.terr[,1:2]


stability1km <- merge(stability1km, coords, by.x = c('eastings', 'northings'), by.y = c('east', 'north'), all.x = FALSE, all.y = TRUE)

## Change recovery time to days rather than units

stability1km$evi.recmed <- stability1km$evi.recmed*8
stability.red <- stability1km[, c(1:3, 5, 9, 11, 13:14)]

stability.red$eastingskm <- stability.red$eastings/1000
stability.red$northingskm <- stability.red$northings/1000

coordskm <- stability.red[, 9:10]



##### Calculating Moran's I

coords.knn <- knearneigh(as.matrix(coordskm), k=8)
coords.nb <- knn2nb(coords.knn)

#create weights
coords.lw <- nb2listw(coords.nb, style='W')
coords.lwB <- nb2listw(coords.nb, style='B')


## Variability

moran(stability.red$evi.var, coords.lw, length(coords.nb), Szero(coords.lw))
#I = 0.533 K = 27.6818

var.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.var)
surface.var <- surf.ls(2, var.df)
variogram(surface.var, 300)
correlogram(surface.var, 300, main = 'Variability')

moran.test(stability.red$evi.var, coords.lw)
#I = 0.533 Exp = -0.000012 Var = 0.0000030 p<0.001

moran.test(stability.red$evi.var, listw = coords.lwB)
#I = 0.533 Exp = -0.000012 Var = 0.0000030


## Magnitude

moran(stability.red$evi.magmed, coords.lw, length(coords.nb), Szero(coords.lw), NAOK = TRUE)
#I =  0.32 K = 16.966

mag.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.magmed)
mag.df <- mag.df[complete.cases(mag.df),]


surface.mag <- surf.ls(2, mag.df)
variogram(surface.mag, 300)
correlogram(surface.mag, 300, main = 'Resistance')

moran.test(stability.red$evi.magmed, coords.lw, na.action = na.omit)
#I = 0.32 Exp = -0.000012 Var = 0.0000030 p<0.001


## Recovery time

moran(stability.red$evi.recmed, coords.lw, length(coords.nb), Szero(coords.lw), NAOK = TRUE)
#I =  0.14 K = 23.0342

rec.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.recmed)
rec.df <- rec.df[complete.cases(rec.df),]

surface.rec <- surf.ls(2, rec.df)
variogram(surface.rec, 300)
correlogram(surface.rec, 300, main = 'Recovery time')

moran.test(stability.red$evi.recmed, coords.lw, na.action = na.omit)
#I = 0.14 Exp = -0.000012 Var = 0.00000301 p<0.001


## Recovery rate

moran(stability.red$evi.recrate.med, coords.lw, length(coords.nb), Szero(coords.lw), NAOK = TRUE)
#I =  0.204 K = 2.489

rate.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.recrate.med)
rate.df <- rate.df[complete.cases(rate.df),]

surface.rate <- surf.ls(2, rate.df)
variogram(surface.rate, 300)
correlogram(surface.rate, 300, main = 'Recovery rate')


moran.test(stability.red$evi.recrate.med, coords.lw, na.action = na.omit)
#I = 0.0204 Exp = -0.000012 Var = 0.00000301 p<0.001


## Longest return

moran(stability.red$evi.long, coords.lw, length(coords.nb), Szero(coords.lw))
#I =  0.29 K = 17.654


long.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.long)


surface.long <- surf.ls(2, long.df)
variogram(surface.long, 300)
correlogram(surface.long, 300, main = 'Longest return')

moran.test(stability.red$evi.long, coords.lw)
#I = 0.29 Exp = -0.000012 Var = 0.00000301 p<0.001


## Slowest return

moran(stability.red$evi.slow, coords.lw, length(coords.nb), Szero(coords.lw))
#I =  0.026 K = 16.443


slow.df <- data.frame(x=stability.red$eastingskm, y = stability.red$northingskm, z = stability.red$evi.slow)

surface.slow <- surf.ls(2, slow.df)
variogram(surface.slow, 300)
correlogram(surface.slow, 300, main = 'Slowest return')

moran.test(stability.red$evi.slow, coords.lw)
#I = 0.026 Exp = -0.000012 Var = 0.00000301 p<0.001



tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\correlogram1km.tiff', height = 14, width =14, unit = 'in',  res = 360)
par(mfrow = c(2,3))
correlogram(surface.var, 300, main = 'Variability', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
correlogram(surface.mag, 300, main = 'Resistance', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
correlogram(surface.rec, 300, main = 'Recovery time', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
correlogram(surface.rate, 300, main = 'Recovery rate', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex =2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
correlogram(surface.long, 300, main = 'Longest return', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
correlogram(surface.slow, 300, main = 'Slowest return', xlim = c(0, 200), xlab = 'Distance', ylab = "Moran's I", cex = 2, cex.main = 2, cex.axis = 2, cex.lab = 2, pch = 19)
dev.off()




