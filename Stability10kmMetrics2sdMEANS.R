#####################################################################################################
### Analyses and graphics for stability metrics at 10 km resolution using 2 sigma events fo means ###
#####################################################################################################

### Hannah White 21.06.2019


library(ggplot2)
library(cowplot)
library(raster)
library(viridis)
library(scales)
library(rasterVis)


## load stability measures

stability10km <- read.csv('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\stability10km2sd.csv', header = TRUE)

## read in location of terrestrial hectads

corine10.terr <- read.csv('J:/Postdoc Grassland Resilience/LandCoverData/corine10.terr.csv', header = TRUE)

coords <- corine10.terr[,2:3]
coords[,1] <- coords[,1] + 5000
coords[,2] <- coords[,2] + 5000

stability10km <- merge(stability10km, coords, by.x = c('eastings', 'northings'), by.y = c('east', 'north'), all.x = FALSE, all.y = TRUE)

## Change recovery time to days rather than units

stability10km$evi.recmed <- stability10km$evi.recmed*8
stability10km$evi.rec <- stability10km$evi.rec*8

stability.red <- stability10km[, c(1:4, 6:7, 13:14)] # variability, mag, rec, recrate, long, slow



## Correlation plot 

library(Hmisc)

rcorr(cbind(stability.red[,3], stability.red[,4], stability.red[,5], stability.red[,6], stability.red[,7], stability.red[,8]),  type = 'spearman') ## 

# pairs plot with histogram
library(psych)

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\correlations10kmhistMEAN.tiff', height = 1000, width = 1000)
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
