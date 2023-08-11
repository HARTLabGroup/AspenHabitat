#Packages Required
#install.packages("tidyverse")
#library(tidyverse)
library(biomod2)
library(raster)
library(rasterVis)
library(ggplot2)
library(gridExtra)
#Heat Load Index
dem.mosaic <- raster(here("C:/Users/paude/Desktop/BioMod2Try/ele30_usgs.tif"))
hli.rast <- hli(dem.mosaic, force.hemisphere = "northern")
writeRaster(hli.rast, here("C:/Users/paude/Desktop/BioMod2Try/-HLI.tif"))

#15 cell neighborhood TPI
dem1 <- raster(here("C:/Users/paude/Desktop/BioMod2Try/ele30_usgs.tif"))
tpi_2<- tpi(dem1, scale=15)
tpi_2
writeRaster(tpi_2, here("C:/Users/paude/Desktop/BioMod2Try/tpi2.tif"), overwrite=TRUE)


#Stacking soil and climate data
setwd("D:/BioMod2Try/0_Data_final/current")
rlist=list.files(pattern="tif$",full.names = TRUE)
s_xvar=stack(rlist)
names(s_xvar)

#Projecting all to WGS_1984 becasue biomod doesn't work in NAD..

#r1=raster("D:/sdm/current_projection/Current_projection.grd")#reference for WGS projection
r1=raster("C:/Users/paude/OneDrive - Colostate/sdm/current_projection/Current_projection.grd")

r1

r2=raster("D:/BioMod2Try/try/clay.tif")#trial if it works

r2

clay1=projectRaster(r2, crs = r1)
clay1


s_xvar1=projectRaster(s_xvar, crs = clay1)#projecting stack downscaled variables along with soil
s_xvar1
Explvar=raster::stack(s_xvar1$ADI,s_xvar1$MSPDD5,s_xvar1$Clay,s_xvar1$pH,s_xvar1$MSP,s_xvar1$om)
#writeRaster(s_xvar1, here("C:/Users/paude/Desktop/BioMod2Try/current.tif")) 

#
Aspen <- read.csv("D:/BioMod2Try/Max_Aspen_Raster/Aspen_Ran2000.csv")
head(Aspen)
summary(Aspen)
#Correlation Maxtix
library(dplyr)#for correlation matrix
presvals_all<-cbind(raster::extract(s_xvar1,Aspen[,c("longitude","lattitude")]),Aspen[,c("longitude","lattitude")])
cor_aspen_all<-presvals_all%>% select(-c("longitude","lattitude"))%>%cor()
cor_aspen_all

write.csv(cor_aspen_all,"D:/BioMod2Try/aspen_cor_finaldata.csv")
##vif
presvals_all<-cbind(raster::extract(Explvar,Aspen[,c("longitude","lattitude")]),Aspen[,c("longitude","lattitude")])
vifcor_1<-presvals_all%>% select(-c("longitude","lattitude"))%>%vifcor()
vifcor_1

class(vifcor_1)
#mat <- as.matrix(vifcor_1)
#print(mat)

#barplot(mat, main = "VIF Values", horiz = TRUE, col = "steelblue")



#Response varaible
AspenData <- read.csv("D:/BioMod2Try/Max_Aspen_Raster/Aspen_Ran2000.csv")
head(AspenData)
myRespName <- "Aspen"
myResp <- as.numeric(AspenData[,myRespName])
myRespXY <- AspenData[,c("longitude","lattitude")]


Explvar=raster::stack(s_xvar1$ADI,s_xvar1$MSPDD5,s_xvar1$Clay,s_xvar1$pH,s_xvar1$MSP,s_xvar1$om)
#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$MSP,s_xvar1$hli,s_xvar1$clay,s_xvar1$ph)

#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$TD)

myExpl=raster::stack(Explvar)


myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName, )

myBiomodData



#plot(myBiomodData)                               
myBiomodOption <- BIOMOD_ModelingOptions()
Aspen_2010 <-
  BIOMOD_Modeling(
    bm.format = myBiomodData,
    models = c("GBM", "RF","GAM"),
    bm.options = myBiomodOption,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 5,
    modeling.id = 'demo1'
  )
Aspen_2010
summary(Aspen_2010)
#get VIF
library(car)
class(Aspen_2010)

vif_values<-vif(Aspen_2010)
vif_vals
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")  
#plot scores
# Get evaluation scores & variables importance

#VarImport <- get_variables_importance(Aspen_models, as.data.frame = TRUE)
VarImport <- get_variables_importance(Aspen_2010)#gives mean of most important variables 

VarImport
apply(VarImport,c(1,2),mean)#mean of imp_variables

bm_PlotEvalMean(Aspen_2010,
                    by="models",
                    metrics=c("ROC","TSS"),
                    main=NULL
)
                
bm_PlotEvalMean(Aspen_2010)

bm_PlotEvalMean(Aspen_2010,
                metrics=c("ROC","TSS")
)

bm_PlotEvalMean(
  Aspen_2010,
  metric.eval = NULL,
  dataset = "calibration",
  group.by = "algo",
  do.plot = TRUE
)
bm_PlotEvalMean(bm.out = Aspen_2010)
#individual models response plots
#Aspen_glm<-BIOMOD_LoadModels(Aspen_2010,models="GLM")
Aspen_rf<-BIOMOD_LoadModels(Aspen_2010,models="RF")
Aspen_gbm<-BIOMOD_LoadModels(Aspen_2010,models="GBM")
Aspen_gam<-BIOMOD_LoadModels(Aspen_2010,models="GAM")

#biomod2::response.plot2( 
glm_eval_strip<-
  bm_PlotResponseCurves(
    bm.out=Aspen_2010,
    models=Aspen_rf,
    Data=get_formal_data(Aspen_2010,'expl.var'),
    show.variables=get_formal_data(Aspen_2010,'expl.var.names'),
    do.bivariate=FALSE,
    fixed.var.metric='median',
    legend=FALSE,
    display_title=FALSE,
    data_species=get_formal_data(Aspen_2010,'resp.var')
  )
glm_eval_strip
##Eval
myEval <- get_evaluations(Aspen_2010, as.data.frame = TRUE)
myEval

myEval$CV.strategy <- "Random"
myEval$CV.strategy[grepl("13", myEval$Model.name)] <- "Full"
myEval$CV.strategy[grepl("11|12", myEval$Model.name)] <- "Stratified"
head(myEval)
p=boxplot(myEval$Testing.data ~ interaction(myEval$Algo, myEval$CV.strategy),
        xlab = "", ylab = "ROC AUC", col = rep(c("brown", "cadetblue"), 3))

p1=bm_PlotVarImpBoxplot(bm.out = Aspen_2010, group.by = c('expl.var', 'algo', 'algo'))

p2=bm_PlotVarImpBoxplot(bm.out = Aspen_2010, group.by = c('expl.var', 'algo', 'dataset'))
p3=bm_PlotVarImpBoxplot(bm.out = Aspen_2010, group.by = c('algo', 'expl.var', 'dataset'))
#impvar= get_variables_importance(Aspen_models,as.data.frame=TRUE)
#plot(impvar)


myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = Aspen_2010,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      prob.mean = TRUE,
                                      prob.median = FALSE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional',
                                      seed.val = 42)
p4=bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
p5=bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'dataset'))
p6=bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'dataset'))

myBiomodProj <- BIOMOD_Projection(bm.mod = Aspen_2010,
                                  proj.name = 'Current',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE,
                                  output.format = '.img')
myBiomodProj

plot(myBiomodProj)

BiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                           proj.name = 'CurrentEM',
                                           new.env = myExpl,
                                           models.chosen = 'all',
                                           metric.binary = 'all',
                                           metric.filter = 'all',
                                           output.format = '.img')

BiomodEMProj
plot(BiomodEMProj)

CurrentProj <- get_predictions(myBiomodProj)
CurrentProj
#"C:/Users/paude/Desktop/BioMod2Try/try/Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
##Species Range Change Required binary stata tp campare, 1 to 1 compasrison between current and furture projection scenario
#Aspen_bin_proj_curr<- (ca="Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
#"C:/Users/paude/Desktop/BioMod2Try/try/Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
Aspen_bin_proj_2010_ca<- (ca="Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
Aspen_bin_proj_2010_wm<- (wm="Aspen/proj_CurrentEM/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")





#a1<-raster("C:/Users/paude/Desktop/BioMod2Try/try/Aspen/proj_CurrentEM/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
#plot(a1)

##Climate Change Scenario SS2_4.5_2010_2040
setwd("D:/BioMod2Try/0_Data_final/LowEmisSce/2040")
#setwd("C:/Users/paude/OneDrive - Colostate/sdm/0_downscaled_variables/Ss4.5_2041_2070")
rlist=list.files(pattern="tif$",full.names = TRUE)
s_xvar=stack(rlist)
names(s_xvar)
#r1=raster("D:/sdm/current_projection/Current_projection.grd")#reference for WGS projection
r1=raster("C:/Users/paude/OneDrive - Colostate/sdm/current_projection/Current_projection.grd")

r1

r2=raster("D:/BioMod2Try/try/clay.tif")#trial if it works

r2

clay1=projectRaster(r2, crs = r1)
clay1


s_xvar1=projectRaster(s_xvar, crs = clay1)#projecting stack downscaled variables along with soil
s_xvar1

#2040
Explvar_2040=raster::stack(s_xvar1$ADI,s_xvar1$MSPDD5,s_xvar1$Clay,s_xvar1$pH,s_xvar1$MSP,s_xvar1$om)
#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$MSP,s_xvar1$hli,s_xvar1$clay,s_xvar1$ph)

#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$TD)

myExpl_2040=raster::stack(Explvar_2040)


myBiomodData_2040 <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl_2040,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName, )
myBiomodData_2040
plot(myBiomodData_2040)                               
myBiomodOption_2040 <- BIOMOD_ModelingOptions()
Aspen_2040 <-
  BIOMOD_Modeling(
    bm.format = myBiomodData_2040,
    models = c("GBM", "RF","GAM"),
    bm.options = myBiomodOption,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 5,
    modeling.id = 'demo1'
  )
Aspen_2040

myBiomodEM_2040 <- BIOMOD_EnsembleModeling(bm.mod = Aspen_2040,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.7),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 2,
                                      prob.mean = TRUE,
                                      prob.median = FALSE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight = TRUE,
                                      prob.mean.weight.decay = 'proportional',
                                      seed.val = 42)

myBiomodProj_2040 <- BIOMOD_Projection(bm.mod = Aspen_2040,
                                  proj.name = '2040',
                                  new.env = myExpl_2040,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE,
                                  output.format = '.img')
myBiomodProj_2040

plot(myBiomodProj_2040)

BiomodEMProj_2040 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_2040,
                                           proj.name = 'myBiomodProj_2040',
                                           new.env = myExpl_2040,
                                           models.chosen = 'all',
                                           metric.binary = 'all',
                                           metric.filter = 'all',
                                           output.format = '.img')


FutureProj_2040 <- get_predictions(myBiomodProj_2040)

b1<-raster("C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2011_2040/Aspen/proj_myBiomodProj/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(b1)
"C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2011_2040/Aspen/proj_myBiomodProj/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"
#"C:/Users/paude/Desktop/BioMod2Try/try/Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
##Species Range Change Required binary stata tp campare, 1 to 1 compasrison between current and furture projection scenario
#Aspen_bin_proj_curr<- (ca="Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

Aspen_bin_proj_2040_wm<- (wm="C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2011_2040/Aspen/proj_myBiomodProj_2040/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
#"C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2011_2040/Aspen/proj_myBiomodProj_2040/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"
Aspen_bin_proj_2040_ca<-(ca="C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2011_2040/Aspen/proj_myBiomodProj_2040/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
Aspen_bin_proj_2040

myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = Aspen_bin_proj_2010_ca, proj.future = Aspen_bin_proj_2040_ca)
SRC<-
  BIOMOD_RangeSize(
    Aspen_bin_proj_2010_wm,
    Aspen_bin_proj_2040_wm
  )
##Climate Change Scenario SS2_4.5_2041_2070
setwd("D:/BioMod2Try/0_Data_final/LowEmisSce/2070")
#setwd("C:/Users/paude/OneDrive - Colostate/sdm/0_downscaled_variables/Ss4.5_2041_2070")
rlist=list.files(pattern="tif$",full.names = TRUE)
s_xvar=stack(rlist)
names(s_xvar)
#r1=raster("D:/sdm/current_projection/Current_projection.grd")#reference for WGS projection
r1=raster("C:/Users/paude/OneDrive - Colostate/sdm/current_projection/Current_projection.grd")

r1

r2=raster("D:/BioMod2Try/try/clay.tif")#trial if it works

r2

clay1=projectRaster(r2, crs = r1)
clay1


s_xvar1=projectRaster(s_xvar, crs = clay1)#projecting stack downscaled variables along with soil
s_xvar1  
####2070
Explvar_2070=raster::stack(s_xvar1$ADI,s_xvar1$MSPDD5,s_xvar1$MSP,s_xvar1$Clay,s_xvar1$pH,s_xvar1$om)
#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$MSP,s_xvar1$hli,s_xvar1$clay,s_xvar1$ph)

#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$TD)

myExpl_2070=raster::stack(Explvar_2070)


myBiomodData_2070 <- BIOMOD_FormatingData(resp.var = myResp,
                                          expl.var = myExpl_2070,
                                          resp.xy = myRespXY,
                                          resp.name = myRespName, )
myBiomodData_2070
plot(myBiomodData_2070)                               
myBiomodOption_2070 <- BIOMOD_ModelingOptions()
Aspen_2070 <-
  BIOMOD_Modeling(
    bm.format = myBiomodData_2070,
    models = c("GBM", "RF","GAM"),
    bm.options = myBiomodOption,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 5,
    modeling.id = 'demo1'
  )
Aspen_2070

myBiomodEM_2070 <- BIOMOD_EnsembleModeling(bm.mod = Aspen_2070,
                                           models.chosen = 'all',
                                           em.by = 'all',
                                           metric.select = c('TSS'),
                                           metric.select.thresh = c(0.7),
                                           metric.eval = c('TSS', 'ROC'),
                                           var.import = 2,
                                           prob.mean = TRUE,
                                           prob.median = FALSE,
                                           prob.cv = FALSE,
                                           prob.ci = FALSE,
                                           prob.ci.alpha = 0.05,
                                           committee.averaging = TRUE,
                                           prob.mean.weight = TRUE,
                                           prob.mean.weight.decay = 'proportional',
                                           seed.val = 42)

myBiomodProj_2070 <- BIOMOD_Projection(bm.mod = Aspen_2070,
                                       proj.name = '2070',
                                       new.env = myExpl_2070,
                                       models.chosen = 'all',
                                       metric.binary = 'all',
                                       metric.filter = 'all',
                                       build.clamping.mask = TRUE,
                                       output.format = '.img')
myBiomodProj_2070

plot(myBiomodProj_2070)

BiomodEMProj_2070 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_2070,
                                                proj.name = 'myBiomodProj_2070',
                                                new.env = myExpl_2070,
                                                models.chosen = 'all',
                                                metric.binary = 'all',
                                                metric.filter = 'all',
                                                output.format = '.img')

FutureProj_2070 <- get_predictions(myBiomodProj_2070)
FutureProj_2070
b1<-raster("C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2041_2070/Aspen/proj_myBiomodProj/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
plot(b1)
"C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ss4.5_2041_2070/Aspen/proj_myBiomodProj/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img"
#"C:/Users/paude/Desktop/BioMod2Try/try/Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
##Species Range Change Required binary stata tp campare, 1 to 1 compasrison between current and furture projection scenario
#Aspen_bin_proj_curr<- (ca="Aspen/proj_CurrentEM/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")

Aspen_bin_proj_2070_wm<- (wm="C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Aspen_2070/proj_myBiomodProj/individual_projections/Aspen_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")

Aspen_bin_proj_2070_ca<-(ca="C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Aspen_2070/proj_myBiomodProj/individual_projections/Aspen_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img"
)

##2100
##Climate Change Scenario SS2_4.5_2071_2100
#C:/Users/paude/Desktop/BioMod2Try/SS245_biomod2/Ssp4.5_2071_2100
setwd("D:/BioMod2Try/0_Data_final/LowEmisSce/2100")
#setwd("C:/Users/paude/OneDrive - Colostate/sdm/0_downscaled_variables/Ss4.5_2041_2070")
rlist=list.files(pattern="tif$",full.names = TRUE)
s_xvar=stack(rlist)
names(s_xvar)
#r1=raster("D:/sdm/current_projection/Current_projection.grd")#reference for WGS projection
r1=raster("C:/Users/paude/OneDrive - Colostate/sdm/current_projection/Current_projection.grd")

r1

r2=raster("D:/BioMod2Try/try/clay.tif")#trial if it works

r2

clay1=projectRaster(r2, crs = r1)
clay1


s_xvar1=projectRaster(s_xvar, crs = clay1)#projecting stack downscaled variables along with soil
s_xvar1  
####2100
Explvar_2100=raster::stack(s_xvar1$ADI,s_xvar1$MSPDD5,s_xvar1$MSP,s_xvar1$Clay,s_xvar1$pH, s_xvar1$om)
#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$MSP,s_xvar1$hli,s_xvar1$clay,s_xvar1$ph)

#Explvar=raster::stack(s_xvar1$ADI1,s_xvar1$GSPDD5,s_xvar1$RH,s_xvar1$TD)

myExpl_2100=raster::stack(Explvar_2100)


myBiomodData_2100 <- BIOMOD_FormatingData(resp.var = myResp,
                                          expl.var = myExpl_2100,
                                          resp.xy = myRespXY,
                                          resp.name = myRespName, )
myBiomodData_2100
plot(myBiomodData_2100)                               
myBiomodOption_2100 <- BIOMOD_ModelingOptions()
Aspen_2100 <-
  BIOMOD_Modeling(
    bm.format = myBiomodData_2100,
    models = c("GBM", "RF","GAM"),
    bm.options = myBiomodOption,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 5,
    modeling.id = 'demo1'
  )
Aspen_2100

myBiomodEM_2100 <- BIOMOD_EnsembleModeling(bm.mod = Aspen_2100,
                                           models.chosen = 'all',
                                           em.by = 'all',
                                           metric.select = c('TSS'),
                                           metric.select.thresh = c(0.7),
                                           metric.eval = c('TSS', 'ROC'),
                                           var.import = 2,
                                           prob.mean = TRUE,
                                           prob.median = FALSE,
                                           prob.cv = FALSE,
                                           prob.ci = FALSE,
                                           prob.ci.alpha = 0.05,
                                           committee.averaging = TRUE,
                                           prob.mean.weight = TRUE,
                                           prob.mean.weight.decay = 'proportional',
                                           seed.val = 42)

myBiomodProj_2100 <- BIOMOD_Projection(bm.mod = Aspen_2100,
                                       proj.name = '2100',
                                       new.env = myExpl_2100,
                                       models.chosen = 'all',
                                       metric.binary = 'all',
                                       metric.filter = 'all',
                                       build.clamping.mask = TRUE,

                                                                              output.format = '.img')
myBiomodProj_2100

plot(myBiomodProj_2100)

BiomodEMProj_2100 <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM_2100,
                                                proj.name = 'myBiomodProj_2100',
                                                new.env = myExpl_2100,
                                                models.chosen = 'all',
                                                metric.binary = 'all',
                                                metric.filter = 'all',
                                                output.format = '.img')

BiomodEMProj_2100
plot(BiomodEMProj_2100)

CurrentProj <- get_predictions(myBiomodProj)
CurrentProj

FutureProj_2100 <- get_predictions(myBiomodProj_2100)

###Species Range Change 2011-2040

SRC_Curr2040 <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj_2040)
SRC_Curr2040

SRC_Curr2040$Compt.By.Models

###Species Range Change 2041-2070

SRC_Curr2070 <- BIOMOD_RangeSize(proj.current = FutureProj_2040, proj.future = FutureProj_2070)
SRC_Curr2070

SRC_Curr2070$Compt.By.Models

###Species Range Change 2011-2100

SRC_Curr2100 <- BIOMOD_RangeSize(proj.current = FutureProj_2070, proj.future = FutureProj_2100)
SRC_Curr2100

SRC_Curr2100$Compt.By.Models



