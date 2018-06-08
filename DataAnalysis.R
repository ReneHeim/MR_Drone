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


# Before we start to analyse the data we set up our working environment.

dir.create('output', FALSE, FALSE)

# install.packages(c("rgdal", 
#                    "raster", 
#                    "roxygen2",
#                    "tictoc", 
#                    "tidyverse", 
#                    "caret",
#                    "e1071", 
#                    "gdata", 
#                    "hsdar", 
#                    "utils", 
#                    "magrittr", 
#                    "rasterVis",  
#                    "rmarkdown"))

library(rgdal) #load installed pkgs
library(raster)
library(tictoc)
library(caret)
library(gdata)
library(hsdar)
library(utils)
library(rasterVis)
library(rmarkdown)
library(roxygen2)
library(magrittr)
library(knitr)
library(reshape2)
library(cowplot)
library(VSURF)


source("R/FUN_raw2speclibhsdar.R")#coverts spec data to hsdar lib
source("R/FUN_drop_cat_var.R")#drops factor and factor level
source("R/FUN_extract_pixel.R")
source("R/FUN_prepggwide2long.R")

# SECTION A --------------------------------------------------------------------

# Loading image brick (1) and 

img <-
  brick("data/20180528_ortho_ground.tif") 

img <- img/65535
  
ndvi <- (img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.3)/(img$X20180528_ortho_ground.5+img$X20180528_ortho_ground.3)
cal <- img$X20180528_ortho_ground.2/img$X20180528_ortho_ground.3
sipi <- (img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.1)/(img$X20180528_ortho_ground.5-img$X20180528_ortho_ground.3)
ari <- (1/img$X20180528_ortho_ground.2)-(1/img$X20180528_ortho_ground.4)
img <- addLayer(img, c(ndvi, cal, sipi, ari))

names(img) <- c("blue", "green", "red", "re", "nir", "alpha", "ndvi", "cal", "sipi", "ari")

# Loading QGIS shape file where sample polygons have been defined (3)
alldata <- shapefile("data/samplepolygons.shp") #3

# Extracting pixel from sample polygons (4) and check data (5)

dfAll <- extract_pixel(img, alldata) #4
head(dfAll) #5

# Use "Treatment" (5) and "ID" (6) which should have been defined in QGIS to 
# construct a new df (7).

Treatment <- alldata$Type #5
ID <- alldata$id #6

dfclass <- cbind(Treatment, ID) #7

# Look up ID and the assigned class (TR, UN, SHD) from the "Type" column (8) and
# then save them in a list (9) while looping through each ID (10).

classvec <- unique(as.character(dfclass[,1])) #8
li <- list() #9

for(i in classvec){
  x <- subset(dfclass, dfclass[,1] == i)
  li[[i]] <- assign(paste0(i,"_","num"), as.numeric(x[,2]))
} #10

# Writes a new df that has written the class labels (11, 12, 14) instead of ID 
# numbers. New df (14) has the "Type" column added (15).

for(i in UN_num){
  dfAll$class[dfAll$class == i] <- "UN" #11
}

for(i in TR_num){
  dfAll$class[dfAll$class == i] <- "TR" #12
}

for(i in SHD_num){
  dfAll$class[dfAll$class == i] <- "SHD" #13
}

dfAll$class <- as.factor(dfAll$class) #14
names(dfAll)[11] <- c("Type") #15

# Divide by 65535 (DN max value) to yield refl between 0 and 1 (16).rename band names according to Micasense RedEdge
# specifications (2). #Band1(blue=475nm), Band2 (green=560nm), Band3 (Red=668nm), 
# Band4 (840nm), Band5 (rededge=717nm), Band6 (alphaband)


dfAll <- dfAll[,c(11,1,2,3,4,5,7,8,9,10)]
names(dfAll) <- c("Type","blue", "green", "red", "re", "nir", "ndvi", "cal", "sipi", "ari")

# Save (17) and reload (18) df containing reflectance values and classes stored 
# in the "Type" column. Also remove alpha (transparency) channel (19) as it is 
# not used for classification. Reorder columns (20) and check if all classes are 
# present (21).

write.csv(dfAll, 'output/2018MyrtleRust_Refl.csv', row.names = FALSE) #17
classif <- read.csv("output/2018MyrtleRust_Refl.csv", check.names = FALSE) #18

# mycor <- cor(classif[,-6]) #test correlation - all bands correlated



# Finally, re-assign (22) and rename (23) column names of the just exportet df 
# to prepare (24) the data for a ggplot2 output (wide to long format). Plot (25)
# spectra to build Figure X.

gplot <- classif[,1:6]

names(gplot) <- c("Type","475", "560", "668", "717", "840")

spectragg <- prep_gg(gplot) #24

a <- ggplot(spectragg, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2) #25

# RF (all classes)--------------------------------------------------------------

# First we set a seed (26) to avoid random number generation in vulnerable 
# processes. Then we partition the extracted pixel data in a training and 
# test subset (27). Finally, we tune the settings for the random forest 
# training process (28, 29).

#set.seed(20180524) #26

inTraining <- createDataPartition(classif$Type, p = .75, list = FALSE)
train <- classif[ inTraining,]
test  <- classif[-inTraining,] #27

rfControl <- trainControl(
  method = "boot",
  number = 100
) #28

rfGrid <- expand.grid(mtry = seq(1, ncol(train)-1, 1)) #29

# Then we initialize a timer (30) to measure the computation time of the 
# training (31) and prediction (32) process. We stop timing (33). Then, we can 
# assign (34) and export (35) a report of our classification. We are also 
# saving (36) the model object to reload (37) and save time when re-running 
# the code.

tic("RF.I") #30

rfFit <- train(Type ~ ., data = train,
               method = "rf",
               importance = TRUE, ntree=500,
               trControl = rfControl, tuneGrid = rfGrid,
               metric = "Accuracy", maximize = TRUE) #31

rfPred <- 
  predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw") #32

toc() #33

Mica.Prediction <- 
  list(fit = rfFit,
       pred = predict.train(rfFit, test[, !names(test) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred, test$Type),
       varImp = varImp(rfFit, scale = TRUE)) #34

sink(file = 'output/I_micadata_report.txt')
Mica.Prediction #35
sink()

saveRDS(Mica.Prediction, 'output/I_micadata_object.rds') #36
Mica.Prediction <- readRDS("output/I_micadata_object.rds") #37

# Feature selection I

fs.1 <- VSURF(classif[,2:10], 
             classif[,1], 
             clusterType = "FORK", 
             ntree = 500,mtry = 4) #warning can be ignored

saveRDS(fs.1, 'output/I_fs_canopy.rds')

fs.norm.1 <- (fs.1$imp.varselect.thres-min(fs.1$imp.varselect.thres))/(max(fs.1$imp.varselect.thres)-min(fs.1$imp.varselect.thres))

out <- classif[,-1]
vi.I <- rbind(names(out[,as.numeric(fs.1$varselect.thres)]),
              round(fs.1$imp.varselect.thres,2), 
              round(fs.norm.1, 2))

row.names(vi.I) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.')


write.csv(vi.I, 'output/vi.I.csv', row.names = TRUE)

# RiskMap ----------------------------------------------------------------------

# Import a brick image (38), similar to the one imported initially (see #1). Any
# ground pixel has been removed (Agisoft Photoscan Pro) to have only the lemon
# myrtle trees available. The image can be used to test the model by predicted 
# the class of each pixel in this images. To nicely display the image, we set 
# background pixel as being tranparent (39). Then, we remove the alpha channel
# (40) which was not used as a predictor variable, divide by the max DN value
# (41) and rename (42) the bands to match the data on which we trained the 
# random forest model. Now we can predict the pixel-classes (43) and plot (44)
# and export the result (45).

imgpred <-  brick("data/20180528_ortho_no_ground.tif")

imgpred <- imgpred/65535

ndvi.p <- (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)/(imgpred$X20180528_ortho_no_ground.5+imgpred$X20180528_ortho_no_ground.3)
cal.p <- imgpred$X20180528_ortho_no_ground.2/imgpred$X20180528_ortho_no_ground.3
sipi.p <- (imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.1)/(imgpred$X20180528_ortho_no_ground.5-imgpred$X20180528_ortho_no_ground.3)
ari.p <- (1/imgpred$X20180528_ortho_no_ground.2)-(1/imgpred$X20180528_ortho_no_ground.4)

imgpred <- addLayer(imgpred, c(ndvi.p, cal.p, sipi.p, ari.p))

names(imgpred) <- c("blue", "green", "red", "re", "nir", "alpha", "ndvi", "cal", "sipi", "ari")

riskpre <- imgpred[[c(1,2,3,4,5,7,8,9,10)]] # remove alpha/transparency channel as it 



riskpred <- predict(riskpre, Mica.Prediction$fit)

plot(riskpred)

currentDate <- Sys.Date()
rstFileName <- paste("output/riskmap",currentDate,".tif",sep="")
writeRaster(riskpred, file=rstFileName, format = "GTiff", bylayer=TRUE, overwrite=TRUE) #45



# SECTION B --------------------------------------------------------------------

# In this section, we code the analysis to compare leaf-scale vs canopy scale 
# classification results. Therefore, we import some leaf spectral reflectance
# data (46) that was recorded on the same plants but in a previous study. Then,
# we drop (47) the class "Healthy" as this was not recorded on the plantation 
# where the aerial imagery was captured. We create a spectral library (48) to 
# make use of the hsdar package and then resample the hyperspectral data to 
# the Micasense RedEdge band specifications (49-52).

hypdata <- 
  read.csv('data/data.wo.out.binned.cut.csv', check.names = FALSE) %>% #46
  drop_class(., .$Type, "Healthy") %>% 
  .[,1:36]#47

speclib2 <- raw2speclib(hypdata) #48

center <-  c(475, 560, 668, 717, 840) #49
fwhm <- c(20, 20, 10, 10, 40) #50
micasense <- as.data.frame(cbind(center, fwhm)) #51

data_mica <- spectralResampling(speclib2, micasense) #52

# Now we would like to use the resampled data (leaf scale) for classification to 
# then compare it with the data that was captured on the canopy level (Micasense
# RedEdge). We extract the data from our spectral library (55) and add a "Type"
# column to assign the classes (56). We also rename the columns according to the 
# Micasense RedEdge bands (57).

micadata <- as.data.frame(data_mica@spectra@spectra_ma) #55
micadata <- cbind('Type'=hypdata$Type, micadata) #56

names(micadata) <- c("Type", "blue", "green", "red", "re", "nir") #57

# We plot the hyperspectral leaf data (53) and also the resampled version (54) 
# for Figure X.

spectraggII <- prep_gg(hypdata)

b <- ggplot(spectraggII, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2) #53

mplot <- micadata
names(mplot) <- c("Type", "475", "560", "668", "717", "840")

spectraggIII <- prep_gg(mplot)

c <- ggplot(spectraggIII, aes(Wavelength, Reflectance, colour = Type)) +
  geom_line(aes(linetype=Type), size = 1)+
  geom_point(aes(shape=Type), size = 2) #54

# Add indices

cal2 <- '560/668'
sipi2 <- '(840-475)/(840-668)'
ari2 <- '(1/560)-(1/717)'

SIPI2 <- vegindex(speclib2, sipi2)
NDVI2 <- vegindex(speclib2, "NDVI")
Calderon2 <- vegindex(speclib2, cal2)
ARI2 <- vegindex(speclib2, ari2)

micadata.ind <- cbind(micadata, NDVI2, SIPI2, Calderon2, ARI2)


# RF (Leaf) --------------------------------------------------------------------

# The seed from our previous classification is still active (#26). We partition 
# the resampled data (58) in a training (59) and test (60) subset. We tune the 
# settings for the random forest training process (61, 62).


inTrainingM <- createDataPartition(micadata.ind$Type, p = .75, list = FALSE) #58
trainM <- micadata.ind[ inTrainingM,] #59
testM  <- micadata.ind[-inTrainingM,] #60

rfControl <- trainControl(
  method = "boot",
  number = 100
) #61

rfGrid <- expand.grid(mtry = seq(1, ncol(trainM)-1, 1)) #62

# Again, we initialize a timer (63) to measure the computation time of the 
# training (64) and prediction (65) process. We stop timing (66). Then, we can 
# assign (67) and export (68) a report of our classification. We are also 
# saving (69) the model object to reload (70) and save time when re-running 
# the code.

tic("RF.II") #63

rfFit.M <- train(Type ~ ., data = trainM,
                 method = "rf",
                 importance = TRUE, ntree=500,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Accuracy", maximize = TRUE) #64

rfPred.M <- #65
  predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw")

toc() #66

mica.leaf.pred <- 
  list(fit = rfFit.M,
       pred = predict.train(rfFit.M, testM[, !names(testM) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred.M, testM$Type),
       varImp = varImp(rfFit.M, scale = TRUE)) #67

sink(file = 'output/II_leaf_report.txt')
mica.leaf.pred #68
sink()

saveRDS(mica.leaf.pred, 'output/II_leaf_object.rds') #69
#mica.leaf.pred <- readRDS("output/II_leaf_object.rds") #70

# Feature selection II

fs.2 <- VSURF(micadata.ind[,2:10], 
              micadata[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

saveRDS(fs.2, 'output/II_fs_leafres.rds')

fs.norm.2 <- (fs.2$imp.varselect.thres-min(fs.2$imp.varselect.thres))/(max(fs.2$imp.varselect.thres)-min(fs.2$imp.varselect.thres))

out2 <- micadata.ind[,-1]
vi.II <- rbind(names(out2[,as.numeric(fs.2$varselect.thres)]),
              round(fs.2$imp.varselect.thres,2), 
              round(fs.norm.2, 2))

row.names(vi.II) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.')


write.csv(vi.II, 'output/vi.II.csv', row.names = TRUE)

# RF (Canopy) ------------------------------------------------------------------

# Now, we have to repeat the classification using the aerial imagery again. This
# time, we drop the class 'SHADOW' (71). This class was not collected for the 
# leaf spectral data and we want to have a fair comparison between leaf-scale 
# and canopy-scale. We can test if the factor levels were dropped (72).

classif2 <- drop_class(classif, classif$Type, 'SHD') #71
unique(classif2$Type) #72


# We can repeat the data partition (73-75), tune the random forest settings 
# (76-77) and then start the training (78) and prediction (79) process. We 
# assign (80) and export (81) a report of our classification. We are also 
# saving (82) the model object to reload (83) if necessary.

inTrainingD <- createDataPartition(classif2$Type, p = .75, list = FALSE) #73
trainD <- classif2[ inTrainingD,] #74
testD  <- classif2[-inTrainingD,] #75

rfControl <- trainControl(
  method = "boot",
  number = 100
) #76

rfGrid <- expand.grid(mtry = seq(1, ncol(trainD)-1, 1)) #77

tic("RF.III")

rfFit.D <- train(Type ~ ., data = trainD,
                 method = "rf",
                 importance = TRUE, ntree=500,
                 trControl = rfControl, tuneGrid = rfGrid,
                 metric = "Accuracy", maximize = TRUE) #78

rfPred.D <- #79
  predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw")

toc()

mica.canopy.pred <- #80
  list(fit = rfFit.D,
       pred = predict.train(rfFit.D, testD[, !names(testD) %in% c("Type")], type = "raw"),
       confusion = confusionMatrix(rfPred.D, testD$Type),
       varImp = varImp(rfFit.D, scale = TRUE))

sink(file = 'output/III_canopy_report.txt')
mica.canopy.pred #81
sink()

saveRDS(mica.canopy.pred, 'output/III_canopy_object.rds') #82
#mica.canopy.pred <- readRDS("output/II_canopy_object.rds") #83

# Feature selection III

fs.3 <- VSURF(classif2[,2:10], 
              classif2[,1], 
              clusterType = "FORK", 
              ntree = 500,mtry = 4)

saveRDS(fs.3, 'output/III_fs_canopy2.rds')

fs.norm.3 <- (fs.3$imp.varselect.thres-min(fs.3$imp.varselect.thres))/(max(fs.3$imp.varselect.thres)-min(fs.3$imp.varselect.thres))

out3 <- classif2[,-1]
vi.III <- rbind(names(out3[,as.numeric(fs.3$varselect.thres)]),
               round(fs.3$imp.varselect.thres,2), 
               round(fs.norm.3, 2))

row.names(vi.III) <- c('Predictor', 'Abs. Imp.', 'Rel. Imp.')


write.csv(vi.III, 'output/vi.III.csv', row.names = TRUE)

# It might be helpful to compare all the relevant spectra.

res <- ggdraw() +
  draw_plot(a, x = 0, y = .66, width = 1, height = .33)+
  draw_plot(b, x = 0, y = .33, width = 1, height = .33) +
  draw_plot(c, x = 0, y = 0, width = 1, height = .33) +
  draw_plot_label(label = c("A", "B", "C"), size = 12,
                  x = c(.12,.12,.12), y = c(0.99,0.66,0.33))

ggsave("output/Figure2.allspectra.png",
       plot = res,
       width = 40,
       height = 20,
       units = "cm",
       dpi = 400
)
#END
