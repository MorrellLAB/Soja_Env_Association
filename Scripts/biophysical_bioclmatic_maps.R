########Plot bioclimatic and biophysical data on maps

library(raster)
library(rgdal)
library(dismo)
library(maptools)
library(maps)

##plot multile maps in same window
library(lattice)
par(mfrow=c(2,2))

#directory where you placed the unziped .bil files for your region of interest
layers_dir<-"directory where you placed your worldclim data"

#get bioclim layers
layers<-list.files(layers_dir,full.names=T,pattern=".bil$")
layers<-lapply(layers,raster)
layers<-stack(layers)
layers_names<-list.files(layers_dir,full.names=F,pattern=".bil$")

#If doing soils, this is the directory where you placed the unziped .tif files for your region of interest
layers_dir<-"directory where soils data live"

#extract soil layer information
layers<-list.files(layers_dir,full.names=T,pattern=".tif$")
layers<-lapply(layers,raster)
layers<-stack(layers)
layers_names<-list.files(layers_dir,full.names=F,pattern=".tif$")


#plot all biophyisical or bioclimatic variables
plot(layers)

# first layer of the RasterStack with ylim being your latitdue and xlim being your longitude
# plot any desired layer by changing the second argument to the appropriate variable, which can be identified in layers_names
# Color scheme can be changed to suit your needs

plot(layers, 1, ylim=c(25, 60), xlim=c(100, 150),col = rev(heat.colors(20))),
col = rev(heat.colors(20)), main="Mean Annual Temperature", xlab="Longitude", ylab="Latitude", add=T)

#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter











