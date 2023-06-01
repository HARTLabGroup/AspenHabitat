library(rgdal)
setwd("D:/BioMod2Try/0_Data_final/postdoc/results/all_clip/EMca/changeinaspen/integer/aspenpointspoly/-1&+1only/withelevpointsextracted/project/samplingtry")
myshapefile <- readOGR(dsn = ".",layer = "proj_ele_High_2040")
myshapefile
class(myshapefile)
my_df<-as.data.frame(myshapefile)
my_df
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),5449), ]
my_sample
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/postdoc/results/all_clip/EMca/changeinaspen/integer/aspenpointspoly/-1&+1only/sampleponts/high_2040_sam5449.csv")

###sample
setwd("D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables")
myshapefile <- readOGR(dsn = ".",layer = "high_2040")
myshapefile
class(myshapefile)
my_df<-as.data.frame(myshapefile)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),5449), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/high_2040.csv")
###high2070
myshapefile1 <- readOGR(dsn = ".",layer = "high_2070")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),5715), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/high_2070.csv")
##high2100
myshapefile1 <- readOGR(dsn = ".",layer = "high_2100")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),6571), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/high_2100.csv")


###low 2040
myshapefile1 <- readOGR(dsn = ".",layer = "low_2040")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),4722), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/low_2040.csv")

##low 2070

myshapefile1 <- readOGR(dsn = ".",layer = "low_2070")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),4637), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/low_2070.csv")

##low 2100
myshapefile1 <- readOGR(dsn = ".",layer = "low_2100")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),5424), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/low_2100.csv")

##Unch_high_2040
myshapefile1 <- readOGR(dsn = ".",layer = "Unch_high_2040")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),86807), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_high2040.csv")

##unch_high_2070
myshapefile1 <- readOGR(dsn = ".",layer = "Unch_high_2070")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),86541), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_high2070.csv")

##unch_high_2100

myshapefile1 <- readOGR(dsn = ".",layer = "Unch_high_2100")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),15181), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_high2100.csv")

##unch_low2040
myshapefile1 <- readOGR(dsn = ".",layer = "Unch_low_2040")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),87534), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_low_2040.csv")

##unch_low2070
myshapefile1 <- readOGR(dsn = ".",layer = "Unch_low_2070")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),87619), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_low_2070.csv")

##unch_low2100
myshapefile1 <- readOGR(dsn = ".",layer = "Unch_low_2100")
myshapefile1
class(myshapefile1)
my_df<-as.data.frame(myshapefile1)
dim(my_df)
my_sample<-my_df[sample(nrow(my_df),86832), ]
write.csv(my_sample,"D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/unch_low_2100.csv")

###Regression
library(MASS)
library(tidyverse)
library(caret)
library(leaps)
install.packages("packagename")

setwd("D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample")
Data <- read.csv("D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/Analysis/allcat.csv")
head(Data)
summary(Data)
p=ggplot(Data,aes(y=Aspen_Category, x=diffMSPDD5))+
  geom_boxplot(position = "dodge")
p
p1=p+stat_summary(fun=mean,geom="points",colour="gray",size=2)
p1
#Binomial reg Specify the predictor variables in the formula
predictors <- c("Elevation","TPI","HLI", "diffADI", "diffMSP", "diffMSPDD5", "clay", "om","pH")  # Replace with your actual predictor variable names

formula <- as.formula(paste("Aspen_binomial~", paste(predictors, collapse = " + ")))

# Fit a multivariate binomial regression model
binomial_model <- glm(formula, family = binomial(), data = Data)

# Print the model summary
summary(binomial_model)
plot(residuals.glm(binomial_model))
abline()
# Perform stepwise model selection using AIC
stepwise_model <- step(binomial_model, direction = "both", scope = formula(binomial_model))
summary(stepwise_model)
stepwise_model$residuals
# Load your dataset (replace "your_dataset.csv" with your actual dataset file)
#data <- read.csv("your_dataset.csv")

models <- regsubsets(Aspen_binomial~Elevation+TPI+HLI+diffADI+diffMSP+diffMSPDD5+clay+om+pH, data = Data, nvmax = 6,
                     method = "seqrep")
summary(models)


##
# Set seed for reproducibility
Apen<-within(Data, Aspen_binomial <- factor(Aspen_binomial, labels = c(0, 1)))
Aspen
ggplot(Aspen, aes(x=Elevation, y=Aspen_binomial, group=1))+geom_point()

install.packages("jtools")
library(jtools)
summ(modelShrub_Pos)

#

predictors <- c("Elevation*diffADI","diffADI", "diffMSP", "pH", "om")  # Replace with your actual predictor variable names

formula <- as.formula(paste("Aspen_binomial~", paste(predictors, collapse = " + ")))

# Fit a multivariate binomial regression model
binomial_model <- glm(formula, family = binomial(), data = Data)
summary(binomial_model)
summ(binomial_model)

###Multinomial Regression
install.packages("nnet")
library(nnet)

Data <- read.csv("D:/BioMod2Try/0_Data_final/plots&results/pointdatawithallvariables/0_sample/Analysis/allcat.csv")
head(Data)
predictors <- c("Elevation*diffADI","diffADI", "diffMSP", "pH", "om")  # Replace with your actual predictor variable names

multinomial <- multinom(Aspen_Category~Elevation+TPI+HLI+diffADI+diffMSP+diffMSPDD5+clay+om+pH, data = Data)
summary(multinom_model)

multinom_model <- multinom(grid_code ~ predictor1 + predictor2, data = your_data)