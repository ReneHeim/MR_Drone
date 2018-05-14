####
# The following code reproduces results of the article:
# 
# Heim, Wright and Oldeland, 2018. Using drones and multispectral imagery to 
#detect myrtle rust on a lemon myrtle plantation. JOURNAL
#
# The code is split in two major sections according to our research questions.
#
# A)  Can we discriminate fungicide treated and untreated lemon myrtle trees 
#     based on multispectral imagery captured by a drone?
#       - Load GIS files to extract raw pixel data for all available spectral 
#         bands and pre-defined classes (TReated, UNtreated, SHAdow)
#       - Raw pixel data will be classified using a random forest classifier and 
#         all spectral bands.
#       - Accuracy metrics will be reported in an error matrix
#
# B)  Is a leaf-scale classification equally accurate as a canopy-scale 
#     classification?
#       - Resample hyperspectral leaf data to multispectral (Micasense) 
#         resolution to be able to compare two classifications based on the same 
#         sensor.
#       - Perform classifcation on this resampled data.
#       - Report accuracies and compare to results found in A) (DN vs refl.)
#
# C)  Are the important spectral regions found on the leaf scale similar to 
#     those found from the air?
#       - Compare relevant spectral regions between A) and B)
####


# Setting up coding environment ------------------------------------------------


dir.create('data', FALSE, FALSE) #creating directories
dir.create('R', FALSE, FALSE)
dir.create('output', FALSE, FALSE)

# install.packages(c("rgdal",
# "raster",
# "tictoc",
# "caret",
# "tidyverse",
# "e1071",
# "randomForest",
# "gdata",
# "hsdar",
# "utils",
# "raster",
# "rasterVis")) #installing required packages (pls uncomment to and run)

library(rgdal) #load installed pkgs
library(raster)
library(tictoc)
library(tidyverse)
library(caret)
library(e1071)
library(randomForest)
library(gdata)
library(hsdar)
library(utils)
library(raster)
library(rasterVis)


source("R/20171224_FUN_raw2speclibhsdar.R")#coverts spec data to hsdar lib
source("R/20170601_FUN_DropCatVar.R")#drops factor and factor level
source("R/FUN_extract_pixel.R")#extracting DN from GTIFF based on a shp file 


# SECTION A --------------------------------------------------------------------


# Load multispectral image

img <- brick("data/FullOrtho_crop.tif") # Change to path that contains .tif file

img@data@names <- c("Blue", "Green", "Red", "NIR", "RedEdge", "Alpha")
plotRGB(img, r = 3, g = 2, b = 1) # test if image is shown correctly, rgb combination arbitrary

# Load data shape

alldata <- shapefile("data/20180427mrdrone_trainpolyforR.shp")

responseCol <- "id"

# Extracting training pixels

dfAll <- extract_pixel(img, alldata)

# Constructing table

Treatment <- alldata$Type
ID <- alldata$id

dfclass <- cbind(Treatment, ID)

# Looks up IDs that are linked to class labels and splits them up

classvec <- unique(as.character(dfclass[,1]))
li <- list()

for(i in classvec){
  x <- subset(dfclass, dfclass[,1] == i)
  li[[i]] <- assign(paste0(i,"_","num"), as.numeric(x[,2]))
  }

# Reassembles a df that has written the class labels instead of ID numbers

for(i in UN_num){
dfAll$class[dfAll$class == i] <- "UN"
}

for(i in TR_num){
  dfAll$class[dfAll$class == i] <- "TR"
}

for(i in SHD_num){
  dfAll$class[dfAll$class == i] <- "SHD"
}

dfAll$class <- as.factor(dfAll$class)

# Rename columns according to Micasense specs

#Band1(blue=475nm), Band2 (green=560nm), Band3 (Red=668nm), Band4 (840nm),
#Band5 (rededge=717nm), Band6 (alphaband)

names(dfAll)[7] <- c("Type")

dfAll[,1:5] <- dfAll[,1:5]/65535#divide by 65535 to yield refl between 0 and 1


# Save and reload table containing raw pixel values and classes

write.csv(dfAll[,1:7], 'output/2018MyrtleRust_Refl.csv', row.names = FALSE)
classif <- read.csv("output/2018MyrtleRust_Refl.csv")
classif <- classif[,c(1,2,3,4,5,7)] # remove alpha band, not required for clas
unique(classif$Type)


# Random Forest Classification -------------------------------------------------


set.seed(20180427)

# Partition Data

inTraining <- createDataPartition(classif$Type, p = .75, list = FALSE)
train <- classif[ inTraining,]
test  <- classif[-inTraining,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
  )

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) 

# RF Model Training

tic("RF") #Start timing

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE)

rfFit #Model output based on training data


# RF Model Testing

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

Mica.Prediction <- 
  list(fit = rfFit,
  pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
  confusion = confusionMatrix(rfPred, test$Type),
  varImp = varImp(rfFit, scale = FALSE))

sink(file = 'output/I_Mica3ClassPred.txt')
Mica.Prediction
sink()

saveRDS(Mica.Prediction, 'output/I_Mica3ClassPred.rds')
Mica.Prediction <- readRDS("output/I_Mica3ClassPred.rds")


# Create RiskMap ---------------------------------------------------------------


imgpred <-  brick("data/20180427_orthophoto_noground_MRDrone.tif") # load lemon 
#myrtle trees without ground to predict wo grass
    riskpre <- subset(imgpred, 1:5) # remove alpha/transparency channel as it 
    #was not used as a predictor var in rf 
        risk <- riskpre/65535 # divide by 65535 to change values to reflectance
            risk@data@names <- 
              c("Blue", "Green", "Red", "NIR", "RedEdge", "Alpha") # rename

            riskpred <- predict(risk, Mica.Prediction$fit)
                plot(riskpred)

                currentDate <- Sys.Date()
                rstFileName <- paste("output/riskmap",currentDate,".tif",sep="")
                writeRaster(riskpred, file=rstFileName, format = "GTiff")
                #NOTE should wrap this file output in function
                
# Visualize RiskMap ------------------------------------------------------------

plotRGB(img, r = 3, g = 2, b = 1)
plot(riskpred, add = T, legend = F)


# SECTION B --------------------------------------------------------------------


# Resample hyperspec to multispec ----------------------------------------------


# Load hyperspectral leaf data

hypdata <- read.csv('data/data.wo.out.binned.cut.csv', check.names = FALSE)
hypdata <- DropClass(hypdata, hypdata$Type, "Healthy")

# Create spectral library to use hsdar pkg

speclib <- raw2speclib(hypdata) #Function needs just numbers as colnames.

# Resample hyperspectral data to Micasense band specifications

center <-  c(475, 560, 668, 717, 840)
fwhm <- c(20, 20, 10, 10, 40)

micasense <- as.data.frame(cbind(center, fwhm))

data_mica <- spectralResampling(speclib, micasense)

plot(data_mica)
plot(speclib)

# Extract reflectance data from micasense spectral library for classification

micadata <- as.data.frame(data_mica@spectra@spectra_ma)
micadata <- cbind('Type'=hypdata$Type, micadata)

# Rename columns

newnamesM <- c("Type", "Blue", "Green", "Red", "RedEdge", "NIR")

names(micadata) <- newnamesM

# Random Forest Classification (MultispecSimulate) -----------------------------

# Partition Data

inTrainingM <- createDataPartition(micadata$Type, p = .75, list = FALSE)
trainM <- micadata[ inTrainingM,]
testM  <- micadata[-inTrainingM,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = seq(1, ncol(trainM)-1, 1)) 

# RF Model Training

tic("RF.M") #Start timing

rfFit.M <- train(Type ~ ., data = trainM,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE)

rfFit.M #Model output based on training data


# RF Model Testing

rfPred.M <- 
  predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

Mica.Resamp.Prediction <- 
  list(fit = rfFit.M,
  pred = predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw"),
  confusion = confusionMatrix(rfPred.M, testM$Type),
  varImp = varImp(rfFit.M, scale = FALSE))

sink(file = 'output/II_ResampLeafPred.txt')
Mica.Resamp.Prediction
sink()

saveRDS(Mica.Resamp.Prediction, 'output/II_ResampLeafPred.rds')


# Random Forest Classification (MultispecDrone) --------------------------------

#Create data (drop shadows)

classif2 <- DropClass(classif, classif$Type, 'SHD')

unique(classif2$Type)

# Partition Data

inTrainingD <- createDataPartition(classif2$Type, p = .75, list = FALSE)
trainD <- classif2[ inTrainingD,]
testD  <- classif2[-inTrainingD,]

# Tune random forest resample process to create variable samples for each tree

rfControl <- trainControl(
  method = "boot",
  number = 100
)

# Approximate mtry (# of vars randomly sampled as candidates at each split)

rfGrid <- expand.grid(mtry = seq(1, ncol(trainD)-1, 1)) 

# RF Model Training

tic("RF.D") #Start timing

rfFit.D <- train(Type ~ ., data = trainD,
                 method = "rf",
                 importance = TRUE, ntree=500,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Accuracy", maximize = TRUE)

rfFit.M #Model output based on training data

# RF Model Testing

rfPred.D <- 
  predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw")

toc() # Stop timing

# Collect results in a list and save

Mica.Sim.Prediction <- 
  list(fit = rfFit.D,
  pred = predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw"),
  confusion = confusionMatrix(rfPred.D, testD$Type),
  varImp = varImp(rfFit.D, scale = FALSE))

sink(file = 'output/III_Drone2UNTRPred.txt')
Mica.Sim.Prediction
sink()

saveRDS(Mica.Sim.Prediction, 'output/III_Drone2UNTRPred.rds')

################################################################################
#NOTES
#
#- Be careful with renaming columns for micasense bands! Check order!! 
#- Compare feature selection between classifications! 
#- How many bootstrap samples are recommended? As many as trees? 
#- Sample size between classifications?