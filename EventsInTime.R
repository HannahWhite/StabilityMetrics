########################################################
### Investigating change in extreme events over time ###
########################################################

### Hannah White 25.04.2019

load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_1km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_1km.RData')


### ggplot barplot
library(ggplot2)


#### All events

## stack values
dd <- evi.dates[which(!is.na(evi.dates))]

tab <- table(cut(dd, 'year'))

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



p <- ggplot(year.df, aes(x = year, y = frequency)) + geom_bar(stat='identity')
p <- p + xlab('Year') + ylab('Frequency')
p <- p + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18))

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\EventsInTime1km.tiff', width = 14, height = 14, units = 'in', res = 360)
p
dev.off()
