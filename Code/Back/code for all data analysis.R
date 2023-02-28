### Soils

```{r, eval=F}
#create new shapefile with wgs_1984 coordinated. Start editing to empty shapefile>draw a polygon masking the SRM, larger than SRM to avoid edge effect>save edit 
#Arctool>sampling>createfishnet>exent newshapefile>#row&#column 10>add field>name "latitude" for Y>add field again name"longitude" for X, note to set in degree decimal> copy paste te record in csv file here as "soildata1.csv"


download.file("https://cran.r-project.org/src/contrib/Archive/XPolaris/")
install.packages(in R studio packages, change from cran to .gz file and upload the downloaded "XPolaris_1.0.2.tar" and then import XPolaris library)
library(curl)
library(sp)
library(sf)
library(tidyverse)
library(terra)
library(spDataLarge)
library(spData)
library(httr)
library(XPolaris)
remove.packages("rlang")
install.packages("rlang")

setwd("C:/Users/aspaudel/OneDrive - Colostate/soildata1")
srm_csv<-read.csv("C:/Users/aspaudel/OneDrive - Colostate/soildata/soillatlong.csv")
View(srm_csv)
xplot(locations = srm_csv[1:2,], localPath = getwd())#try XY coordinate locates in map
xplot(locations = srm_csv, localPath = getwd()) #XY of study area
df.loc<-srm_csv
df.loc
xplot(locations = df.loc,localPath = getwd())
#downloading as Polaris output, soil data in various depth in unit of cm
df.out <- ximages(locations = srm_csv,
                  variables = c('ph','om','clay','theta_s'), 
                  statistics = c('mean'), # Only the mean will be 
```


###BioMod2

library(raster)
library(ggplot2)
library(sf)
library(dplyr)
##Presence data
DataAspen <- read.csv("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/aspen_thinned_1t.csv")
head(DataAspen)
##Env data projection of wgs1984 (uniform to aspen presence data)#current predictor
try1<-raster("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/xy_wgs1984.tif")
try2<-raster::stack("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/Joshenvdata/predictors_current.tif")
Env_1980<-projectRaster(try2,crs=crs(try1))
env<-stack(Env_1980)

##Env data projection of wgs1984 #2040
try3<-raster::stack("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/Joshenvdata/predictors_370_2040.tif")
Env_2040<-projectRaster(try3,crs=crs(try1))
env_2040<-stack(Env_2040)
##Env data projection of wgs1984 #2070
try4<-raster::stack("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/Joshenvdata/predictors_370_2070.tif")
Env_2070<-projectRaster(try4,crs=crs(try1))
env_2070<-stack(Env_2070)
##Env data projection of wgs1984 #2100
try5<-raster::stack("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/Joshenvdata/predictors_370_2100.tif")
Env_2100<-projectRaster(try5,crs=crs(try1))
env_2100<-stack(Env_2100)
##
DataAspen <- read.csv("C:/Users/paude/OneDrive - Colostate/Postdoc_CSU/Data/aspen_thinned_1t.csv")
head(DataAspen)
myRespName <- 'Aspen'
myResp <- as.numeric(DataAspen[,myRespName])
myRespXY <- DataAspen[,c("longitude","lattitude")]

library(raster)

myExpl<-env
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'random')

myBiomodOption <- BIOMOD_ModelingOptions()
