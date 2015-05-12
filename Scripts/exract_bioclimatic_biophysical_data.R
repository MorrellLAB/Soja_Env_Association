#Packages to extract the data
library(raster)
library(rgdal)
require(raster)

#download the soil property layers from http://soilgrids1km.isric.org/index.html
#you only need the BLD (bulk density), CEC (cation exchange), CLYPPT (percent clay), ORCDRC (soil organic matter)
#PHIHOX (soil PH), SLTPPT (percent clay), SNDPPT (percent sand)
#you get 6 depths 0-5cm, 5-15cm, 15-30cm, 30-60cm, 60-100cm, 100-200cm
# you can download upper, lower or mean estimates for each of these variables
#This can be used on worldclim data as well the key is to make sure you have the correct file extension 
#file extentions for worldclim are .bil 

#directory where you placed the unziped .tif files for your region of interest
layers_dir<-"directory where your data live"

#file that contains the latitute and longitude of your study sites
points<-read.table("file that contains latitute and longitude of your sampling locations")

points2<-read.table("file that contains latitute and longitude of your sampling locations")

#extract only longitude and latitude columns from your data file
coordinates(points)<-c("longitude","latitude")

#make your program know your points are spatially explicit
points<-SpatialPoints(points)

#extract soil layer information

layers<-list.files(layers_dir,full.names=T,pattern=".tif$")
layers<-lapply(layers,raster)
layers<-stack(layers)
layers_names<-list.files(layers_dir,full.names=F,pattern=".tif$")

#extract values from the layers
values<-extract(layers,points)

#clean up values
final<-cbind(points2,values)
names<-colnames(final)
names<-sub(pattern ="_M_02_apr_2014",replacement = "",x =names )
colnames(final)<-names

#put soil values in a csv file
write.csv(final,"file name of your final csv",row.names=F)

